---
layout: post
title: "Eight Places I Got Stuck Setting Up SSH from a Mac to a Home Windows PC"
date: 2026-12-23 20:00:00 +0900
lang: en
lang-exclusive: ["en"]
permalink: /posts/mac-windows-ssh-setup/
page_id: mac-windows-ssh-setup
categories: [Notes]
tags: [ssh, windows, openssh, ddns, port-forwarding]
description: "Eight traps hit while building an SSH path from a Mac to a home Windows desktop, and how six of them came down to looking at the wrong machine or the wrong file."
---

I set up a working SSH path from my Mac to the Windows desktop at home. I expected
"install OpenSSH, open a port, done." In practice I got stuck again and again, in places
I had not anticipated. This write-up weighs **where things broke and why** more heavily than
the configuration that finally worked.

For the basic job of putting an OpenSSH server on Windows, I followed
[a Korean walkthrough](https://blog.naver.com/twilight_teatime/223668641470) and
[Microsoft's install documentation](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_install_firstuse).
What follows is everything after that.

> Every IP address, domain and port below is an example value. Substitute your own.

---

## The final setup

```
[Mac]  ──internet──▶  router (public IP)  ──port forward──▶  [Windows PC]
                      external port 1927       192.168.0.59:4131
                              │
                       DDNS (No-IP) keeps
                        myhost.ddns.net
```

- Windows: OpenSSH Server on port `4131` (moved off the default 22)
- Router: external `1927` → internal `192.168.0.59:4131`
- DDNS: a dynamic public IP pinned to `myhost.ddns.net`
- Authentication: ed25519 public key
- Clients: the macOS terminal and VS Code Remote-SSH

---

## The eight traps

The conclusion first. There were eight.

| # | Symptom | Actual cause |
|---|---|---|
| 1 | Password prompt persists after installing a public key | Administrator accounts do not read `authorized_keys` |
| 2 | Key registered, authentication still fails | **Pasted the placeholder verbatim** from the instructions |
| 3 | `pbcopy: command not found` | Ran a macOS command inside the SSH session |
| 4 | Every earlier setting was void | The router had changed and the whole subnet with it |
| 5 | Wondering whether to pick DynDNS | Its free plan was retired long ago |
| 6 | Router DDNS authentication fails | The DDNS Key method wants a different value in the domain field |
| 7 | Forwarding configured, still timing out | Tested from inside the house and hit NAT loopback |
| 8 | Terminal works, only VS Code fails | Pasted the whole `ssh` command into the host field |

---

## Trap 1. Administrator accounts do not read `authorized_keys`

The first and longest detour.

I put the public key in `C:\Users\user\.ssh\authorized_keys` and still got a password prompt.
Permissions were right, the key format was fine. The answer was at the bottom of `sshd_config`.

```
Match Group administrators
       AuthorizedKeysFile __PROGRAMDATA__/ssh/administrators_authorized_keys
```

OpenSSH for Windows reads **a separate shared file for accounts in the Administrators group**.
Mine was an administrator account, so the `authorized_keys` in my home directory was never
consulted at all.

The real location is this:

```
C:\ProgramData\ssh\administrators_authorized_keys
```

To check whether your account is an administrator:

```powershell
Get-LocalGroupMember -Group Administrators
```

Permissions matter too. This file must be accessible to **Administrators and SYSTEM only**,
or sshd silently ignores it. No error is emitted.

```powershell
icacls C:\ProgramData\ssh\administrators_authorized_keys /inheritance:r /grant "Administrators:F" /grant "SYSTEM:F"
Restart-Service sshd
```

This behavior is spelled out in Microsoft's
[OpenSSH key management documentation](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_keymanagement).
Reading it first would have saved the whole detour.

> **Lesson:** when public-key authentication fails on Windows, suspect the file location first.
> Coming from Linux only, this is not a behavior you would predict.

---

## Trap 2. I pasted the placeholder verbatim

The instructions contained this command:

```powershell
Add-Content -Path C:\ProgramData\ssh\administrators_authorized_keys -Value 'ssh-ed25519 AAAA...pasted-content... eddy@mac' -Encoding ascii
```

The real key belongs where `AAAA...pasted-content...` sits. I **ran the line as written**.
Authentication failed, of course, and a junk line was left in the file.

The interesting part is that **a file-size change alone** caught it.

```
before: 749 bytes
after:  787 bytes   → +38
```

An ed25519 public key line normally runs around 90 bytes. Thirty-eight is far too small.
And with `-Encoding ascii`, each non-ASCII character collapses to a single `?`. Count the
placeholder that way and you get exactly 38 bytes: 12 for `ssh-ed25519` and a space,
15 for the placeholder token, 1 space, 8 for `eddy@mac`, and 2 for the `CRLF`.

The cleanup:

```powershell
$p = 'C:\ProgramData\ssh\administrators_authorized_keys'
(Get-Content $p) | Where-Object { $_ -notmatch '\.\.\.' } | Set-Content $p -Encoding ascii
```

> **Lesson:** after registering a public key, **check the file size.**
> Around +90 bytes for ed25519, around +750 for RSA 4096.
> sshd skips a malformed key line silently, so no error will tell you.

---

## Trap 3. I ran macOS commands inside the SSH session

The order for creating a key on the Mac is this:

```bash
ssh-keygen -t ed25519 -C "eddy@mac"
pbcopy < ~/.ssh/id_ed25519.pub
ssh -p 4131 eddy@192.168.0.59
```

I had **already connected over SSH** before running them. The results:

- `ssh-keygen` created a key on Windows, not on the Mac
- `pbcopy` is macOS-only, so Windows reported an unknown command

After an SSH connection the terminal window still belongs to the Mac, but **every keystroke
goes to the remote machine.** I knew that, and the identical window kept fooling me anyway.

Now, whenever in doubt, I type this:

```bash
uname -s
```

- `Darwin` → the Mac
- an error or anything else → the remote Windows machine

The real fix was structural: **keep two terminal tabs, one local-only and one SSH-only.**

---

## Trap 4. Replace the router and every old note becomes garbage

An old note of mine said:

```
203.0.113.17:1927 -> 192.168.55.27:4131
```

Checking the current values, all of them had changed.

```
public IP    203.0.113.17    →  203.0.113.88
internal IP  192.168.55.27   →  192.168.0.59
gateway      192.168.55.1    →  192.168.0.1
```

When a move or a new router **shifts the whole subnet**, the forwarding rules simply do not
exist on the new device. Re-read the current values before starting, or the work is wasted.

```powershell
Get-NetIPAddress -AddressFamily IPv4 | Where-Object { $_.PrefixOrigin -eq 'Dhcp' }
Get-NetRoute -DestinationPrefix '0.0.0.0/0' | Select-Object NextHop
(Invoke-RestMethod 'https://api.ipify.org?format=json').ip
```

### Do this first: check for CGNAT

Before spending any time on DDNS and port forwarding, settle this.

Compare the WAN IP on the router's admin page against the public IP you just queried.

| Router WAN IP | Verdict |
|---|---|
| Same as the public IP | Fine. Port forwarding will work |
| `100.64.x.x` – `100.127.x.x` | CGNAT. Port forwarding is impossible |
| `10.x.x.x`, `192.168.x.x`, `172.16–31.x.x` | Double NAT. Another router sits upstream |

In the latter two cases neither port forwarding nor DDNS will work. Either ask the ISP for a
public IP, or **switch to an overlay VPN such as Tailscale.** Honestly that route is easier:
no router configuration, and no port exposed to the internet. The only argument for port
forwarding is not wanting to install extra software.

My WAN IP matched the public IP, so I carried on.

---

## Trap 5. DynDNS is no longer an option

The router's DDNS page offered exactly two providers: `NO-IP` and `DynDNS`.

But **DynDNS retired its free hostname service in May 2014.** Per Dyn's own notice, free
accounts were disabled at that point and had to move to the paid Remote Access product
(then called DynDNS Pro). Oracle runs that paid service today. The entry surviving in router
firmware is just a stale list.

So I chose **No-IP**, whose free plan allows up to three hostnames.

### Creating a No-IP hostname

`My No-IP → Dynamic DNS → Hostnames → Create Hostname`

| Field | Value |
|---|---|
| Type | `A` |
| Host | any name (must be globally unique) |
| Domain | a free domain such as `ddns.net` |
| IPv4 | the current public IP, filled in automatically |
| TTL | Auto |
| Wildcard | off (paid feature) |
| **Enable Dynamic DNS** | **checked** |

The last field is the one that matters. Leave it unchecked and you get a **static A record**,
which the router cannot update when the IP changes — removing the entire reason to use DDNS.

### The free plan's catch

Free No-IP requires **clicking a confirmation email link every 30 days**, or the hostname is
deleted. If something that worked for months suddenly stops, this is almost always why.

### Choosing a domain

Widely used free DDNS domains such as `ddns.net` and `hopto.org` have been reported blocked by
corporate firewalls and public Wi-Fi DNS filters. I have not measured this across networks
myself, so I would not generalize it — but if you expect to connect from an office, it is worth
knowing the name itself can be what fails. No-IP offers several domains even on the free plan,
so recreating the hostname under a different one is an easy fallback.

---

## Trap 6. The DDNS Key method wants a different value in the domain field

Checking `Enable Dynamic DNS` leads to a screen that issues a **DDNS Key**. This is the safer
route, since the router never has to store your No-IP account password.

Here is where I got confused. The issuing screen says:

> Use Hostname for Integrated DDNS Devices: **all.ddnskey.com**

So the values to enter on the router are these:

| Router field | Value |
|---|---|
| Service provider | `NO-IP` |
| Username | the DDNS Key username (`xxxxxxx@ddnskey.com` form) |
| Password | the DDNS Key password |
| Domain name | **`all.ddnskey.com`** ← not your own hostname |
| Bind to WAN IP | **off** |

`all.ddnskey.com` is a special address meaning "update every hostname attached to this key."
The router signals that address and No-IP updates `myhost.ddns.net` on its own. No-IP's
[router DDNS guide](https://www.noip.com/support/knowledgebase/how-to-configure-ddns-in-router)
gives the same instruction for the domain field.

**The address you connect to is still `myhost.ddns.net`.** `all.ddnskey.com` goes only in the
router's configuration field. This was the single most confusing part.

> In practice, entering the hostname (`myhost.ddns.net`) directly also worked.
> If the DDNS Key is assigned to that hostname, either value is accepted.
> Switch to `all.ddnskey.com` if authentication fails.

**The password is shown exactly once, at issue time.** Open a note first, then press create.
Reissuing it means redoing the router configuration from scratch.

### Leave "bind to WAN IP" off

Turn it on and the router reports **its own WAN interface IP** to No-IP. Under CGNAT that
registers a private address in DNS and **breaks a record that was working fine.**

Left off, No-IP auto-detects the source IP of the incoming packet — the real public address.
For an ordinary home router, off is correct.

---

## Trap 7. Testing from inside the house fails

I configured the forwarding rule. The menu path below is for a TP-Link router.

```
Advanced → NAT Forwarding → Virtual Servers → Add
  external port  1927
  internal IP    192.168.0.59
  internal port  4131
  protocol       TCP
```

> **Port Triggering** in the same menu is easy to confuse with this.
> **Virtual Servers** is the right one.

Testing from the Mac on the house Wi-Fi gave a timeout. I combed through the rules assuming
they were wrong; the cause was **NAT loopback**, also called hairpin NAT.

That is the path where traffic leaves the internal network toward its own public IP and comes
back in, and **many routers do not support it.** Perfect forwarding can still fail from inside.

**Test from an external network — a phone hotspot, with the Mac's Wi-Fi turned off.**
To check only whether the port is actually open, an external port checker such as
[canyouseeme.org](https://canyouseeme.org/) does the job.

Pin the internal IP while you are in there, too. If DHCP reassigns it, forwarding breaks silently.

```
Advanced → Network → DHCP Server → Address Reservation
```

---

## Trap 8. The terminal works, only VS Code fails

The terminal connected fine, but VS Code Remote-SSH produced this:

```
WARN: $PLATFORM is undefined in installation script output. Errors may be dropped.
Failed to parse remote port from server output
Exec server for ssh-remote+7b22686f73744e616d65223a226d79686f73742e64646e732e6e6574222c...
```

The message looks unhelpful, but **the trailing hex is just JSON encoded as hex.** Decoded:

```json
{"hostName":"myhost.ddns.net","user":"ssh -p 1927 eddy"}
```

There is the culprit. I had pasted `ssh -p 1927 eddy@myhost.ddns.net` whole into the host
field, so **`ssh -p 1927 eddy` became the username.** The port fell back to 22, which could
never work.

To decode:

```bash
echo '7b22686f73744e616d65...' | xxd -r -p
```

> **Lesson:** always decode the hex in a Remote-SSH error. The cause is written inside it.

### The fix: write it into `~/.ssh/config`

Typing into the host field invites this mistake repeatedly. Write it once in the **Mac's**
`~/.ssh/config` and it is settled.

```
Host winpc
    HostName 192.168.0.59
    User eddy
    Port 4131
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 30

Host winpc-ext
    HostName myhost.ddns.net
    User eddy
    Port 1927
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 30
```

```bash
chmod 600 ~/.ssh/config
```

The terminal gets shorter,

```bash
ssh winpc        # at home
ssh winpc-ext    # away
```

and in VS Code, `Cmd+Shift+P → Remote-SSH: Connect to Host...` now **lists `winpc` and
`winpc-ext` directly.** Nothing left to type.

> You can open the file from VS Code via
> `Cmd+Shift+P → Remote-SSH: Open SSH Configuration File...`.
> The first connection takes a minute or two while VS Code installs its server component
> on the remote machine.

---

## Closing checklist

Done in this order, most of the above is avoidable.

1. **Check for CGNAT** — compare the router WAN IP against the public IP. If they differ, stop here and switch to Tailscale
2. Install OpenSSH Server on Windows; if you change the port, **add the firewall rule for that port too**
3. On an administrator account, the key goes in `C:\ProgramData\ssh\administrators_authorized_keys`
4. After registering the key, **check the size increase** (around +90 bytes for ed25519)
5. Get a connection working **inside the same network** first
6. Create the No-IP hostname with `Enable Dynamic DNS` checked, then configure router DDNS with WAN IP binding off
7. Port forwarding plus a DHCP address reservation
8. **Switch to a hotspot** and test from outside
9. Write the Mac's `~/.ssh/config`
10. In VS Code Remote-SSH, pick the host from the list

## Diagnostic commands

On the Windows side:

```powershell
Get-Service sshd                                             # service state
Get-NetTCPConnection -State Listen -LocalPort 4131           # listening check
Get-Item C:\ProgramData\ssh\administrators_authorized_keys   # key file size
Resolve-DnsName myhost.ddns.net -Type A -DnsOnly             # DDNS propagation
(Invoke-RestMethod 'https://api.ipify.org?format=json').ip   # current public IP
```

On the Mac:

```bash
ssh -v -p 1927 eddy@myhost.ddns.net    # trace where it stops
dig +short myhost.ddns.net             # DNS lookup
uname -s                               # local or remote
```

Reading the cause out of `ssh -v`:

| Message | Cause |
|---|---|
| `Operation timed out` | Missing port forward, or NAT loopback |
| `Connection refused` | Forwarding works, the internal port is wrong |
| `Could not resolve hostname` | DDNS not propagated, or a typo |
| `Permission denied (publickey,password)` | **The network works.** Only key registration is left |

That last line means you are halfway there: the network path is open and only authentication
remains.

---

## Securing it afterwards

Once the connection works, a few closing items.

**1. Disable password authentication**

Only **after** confirming that public-key access works, and with one connected session left
open, make the change. Copying `sshd_config` before editing it is the safer habit.

```
# C:\ProgramData\ssh\sshd_config
PasswordAuthentication no
```

The moment a port is open to the internet, scans arrive from everywhere. Keep password
authentication closed.

**2. Do not leave credentials in plaintext**

It is easy to jot a key passphrase or a DDNS Key password into the notes file you made while
setting this up. If the PC is compromised those go with it. Move them to a password manager.

**3. Change the default shell (optional)**

OpenSSH on Windows defaults to `cmd.exe`. If you prefer PowerShell:

```powershell
New-ItemProperty -Path "HKLM:\SOFTWARE\OpenSSH" -Name DefaultShell `
  -Value "C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe" -PropertyType String -Force
```

That said, VS Code Remote-SSH behaves most reliably under `cmd.exe`. If Remote-SSH is your
main client, leaving it alone is the better choice.

---

## Looking back

Of the eight traps, **six were not "the setting was wrong" but "I was looking in the wrong place."**

- Key in the wrong file (trap 1)
- Wrong content in it (trap 2)
- Commands typed on the wrong machine (trap 3)
- Old values from a different network (trap 4)
- Hostname in the wrong field (trap 6)
- Testing from the wrong location (trap 7)

SSH configuration is not hard. The real difficulty is that the setup makes it easy to lose
track of **which file on which machine you are currently looking at.**

Two habits helped most:

- **In doubt, run `uname -s`** — local or remote?
- **After changing something, check the size and the lookup result** — plenty of settings are ignored silently, with no error

---

## References

- [Install OpenSSH for Windows](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_install_firstuse) — Microsoft Learn
- [OpenSSH key management](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_keymanagement) — Microsoft Learn; the source for `administrators_authorized_keys` behavior
- [Setting up an SSH server on Windows](https://blog.naver.com/twilight_teatime/223668641470) — the basic install walkthrough (Korean)
- [Changing the OpenSSH server port](https://www.lainyzine.com/ko/article/how-to-change-port-number-of-openssh-server-on-windosw-10/) — lainyzine (Korean)
- [Free user discontinuation FAQ](https://help.dyn.com/free-user-discontinuation-faqs/) — Dyn; the source for when the free plan ended
- [How to configure DDNS in a router](https://www.noip.com/support/knowledgebase/how-to-configure-ddns-in-router) — No-IP
- [How to set up and use a DDNS Key](https://www.noip.com/support/knowledgebase/how-to-setup-and-use-a-ddns-key) — No-IP
- [DDNS on a TP-Link router](https://www.noip.com/support/knowledgebase/setting-ddns-tp-link-router) — No-IP
- [Remote Development using SSH](https://code.visualstudio.com/docs/remote/ssh) — VS Code documentation
- [Open Port Check Tool](https://canyouseeme.org/) — verify a port is open from outside
