---
layout: post
title: "Six Places I Got Stuck Driving a MacBook and a Windows Desktop from One Keyboard"
date: 2026-12-25 20:00:00 +0900
lang: en
lang-exclusive: ["en"]
permalink: /posts/synergy-mac-windows-setup/
page_id: synergy-mac-windows-setup
categories: [Notes]
tags: [synergy, autohotkey, macos, windows, clamshell]
description: "Six walls hit while driving a clamshell MacBook and a Windows desktop from one keyboard and mouse with automatic monitor switching. Three are still unsolved."
---

I built a setup where one keyboard and mouse drive both a MacBook (clamshell, on a Satechi
cube dock) and a Windows desktop, with two monitors divided between them automatically
depending on what I am doing. It looked simple. There were about six unexpected walls, and
this is the record of them.

> This is a November 2026 snapshot. Three of the six — monitor flicker, clamshell display
> wake, and automatic reversion — are still unsolved, and each section states its current status.

## The setup

| Item | Detail |
| --- | --- |
| Mac | MacBook Pro (M4 Pro), clamshell mode, on a Satechi cube dock |
| Windows | Desktop (RTX 4080), **Synergy primary (server)** |
| Monitor 1 | ASUS VA24DQLB — Windows only |
| Monitor 2 | LG UltraGear — shared between Windows and the Mac via input switching |
| Software | Synergy 3 Lite v3.6.3, AutoHotkey v2, ControlMyMonitor, KeepingYouAwake |

### Physical wiring

```
[Windows desktop] ──(DP/HDMI)──> [VA24DQLB]  … Windows only
        │
        └────────(HDMI)────────> [LG UltraGear] ─── input 3
        │
        └──────(direct ethernet)─> [Satechi cube dock] ── [MacBook Pro]
                                          │
                                          └──(HDMI)──> [LG UltraGear] ─── input 1
```

The essentials:

- **The two machines are wired directly to each other** — a 1:1 link with no router in
  between, using the ethernet port on the Satechi dock for the Mac side.
- **The desktop holds both monitors; the MacBook holds only the LG.** So the LG has cables
  from both computers at once, and input switching decides which one owns it.
- When the Mac is idle, the LG serves as an extended display for Windows. When the Mac is in
  use, the LG switches to the Mac and Windows falls back to a single monitor.

> The direct ethernet link has no DHCP, so both ends take `169.254.x.x` (APIPA) addresses.
> Synergy connects over those; the SSH control added later uses the Wi-Fi addresses that go
> through the router. APIPA addresses can change on reboot, so pinning static IPs is safer.

---

## Wall 1. Synergy refuses to connect

### Symptom

Adding the Windows IP from the Mac produced **"Serial key mismatch."**

### What I tried, and what it told me

- I suspected the network, but the very fact that this error appeared was **proof that the
  other computer had been reached and had answered.** A failed connection would have said
  "host unreachable."
- The cause was licensing. The Mac was on **Synergy 3 Max** and Windows was using a different key.

### Fix

Unify both ends on **Synergy 3 Lite**:
About → Update serial key → sign in with the same Symless account.

### Then a second failure

With the keys matched, the error became **"Failed."**
The IP I had entered was `192.254.214.194` — a typo for `169` in the first octet.

> **Lesson**: a changed error message means the previous step is fixed.
> Stop suspecting the old cause and read the new message as new.

### Result

Corrected to `169.254.214.194`, connected, mouse and keyboard shared.

---

## Aside: the Synergy license tiers

A question that came up while matching serial keys. Both ends need the same key anyway, and I
do not use Synergy's own hotkeys, so is there any reason to go above the base tier?

Per Symless's documentation, the personal tiers — Lite, Plus and Max — support 3, 5 and 15
computers respectively as one-time purchases, and **TLS/SSL network encryption and hotkeys are
in Plus and Max only.** As of November 2026 the prices are $15, $19 and $29. Copy and paste
between machines is in every tier.

So the **Lite tier I use is the minimum: no encryption, no hotkeys.** For two computers with
no need for Synergy's hotkeys, there is no reason to move up.

The free open-source alternative is **Deskflow**. Interestingly, Deskflow is Synergy's official
upstream project, and Synergy sponsors it with both code and funding. Anyone starting fresh
would do well to try it before paying.

---

## Aside: when clipboard sharing does not work

**Symptom**: mouse and keyboard sharing work fine, but text copied on Windows will not paste
on the Mac.

**What to check**

1. Synergy → Settings → Copy & paste toggle (on by default)
2. macOS System Settings → Privacy & Security → **Input Monitoring** → add Synergy

**Cause**: the license change reinstalled the app, which reset its permissions. The keystroke
permission (Accessibility) survived, so mouse and keyboard kept working perfectly — leaving
only the clipboard broken, which made the cause hard to guess.

---

## Wall 2. The right Alt key is invisible to the Mac

### Background: a Windows-first user picking up a Mac

For someone who has always used a Mac this is a non-issue. Mac keyboards have a dedicated
input-source key, and Caps Lock switches too.

The problem is **a case like mine: a Windows-first user who added a Mac.** On a Korean
keyboard the right Alt acts as the Hangul/English key (VK_HANGUL) under Windows, and macOS has
no concept that corresponds to it.

I already had a workaround for the Mac on its own, using **Karabiner-Elements**. I set the
macOS input-switch shortcut to **F18** and mapped **right Command → F18** in Karabiner. On an
external keyboard the right Command sits exactly where the right Alt is, so physically it came
out as "right Alt switches input source."

### Symptom

With Windows as the primary and Synergy connected, **that Karabiner mapping did nothing at all.**

### Why it fails

Two reasons stacked on top of each other.

1. **Karabiner works by grabbing physical HID devices** and intercepting their input. The
   keystrokes Synergy delivers to the Mac do not come from a physical device; they are
   **events synthesized and injected by software**, so they pass straight through the
   Karabiner layer. The rule never fires at all.
   Leave Karabiner-EventViewer open and press the key: nothing is logged.

2. **Synergy's modifier mapping cannot distinguish left from right.** Write `alt` in the
   config and it applies to both Alt keys, so "right Alt only" is structurally impossible at
   the config level.

### Fix: move the translation to the Windows side

I changed the macOS input-switch shortcut **from F18 to F16** — F18 was already taken by
Karabiner and I wanted to avoid the collision — and switched to **firing F16 directly from
Windows with AutoHotkey.**

Synergy passes that F16 through to the Mac, and macOS system shortcuts respond to synthesized
events perfectly well. Karabiner is not in the path at all.

```autohotkey
~vk15::Send "{Blind}{F16}"
```

- `vk15` = the virtual key code of the Hangul/English key
- `~` = let the original key through as well, so **Windows input switching still works** and
  F16 is emitted in addition
- With the cursor on Windows, F16 does nothing, so there is no side effect

I left the Mac-only Karabiner mapping (right Command → F18) in place, since it is still valid
when using the Mac's own keyboard.

### Result

One key switches input source on both operating systems.

---

## Wall 3. The LG monitor flickers or cuts out intermittently

### What I tried

I suspected Synergy's `screenSaverSync = false`; it was unrelated. Synergy has nothing to do
with video output.

### A partial cause

The LG monitor's **deep-sleep power mode.** To cut standby draw it shuts down the
signal-detection circuit as well, so when the dock powered back up the monitor never noticed
the signal. Turning that option off reduced the "dock is on but nothing appears" symptom.

For reference, **Smart Energy Saving** in the same menu adjusts brightness to screen content
and is unrelated to this.

### What remains

The momentary flicker during use is **still there.** Three candidates:

- **Insufficient dock power**: on a USB-C dock running on bus power alone, with no PD charger,
  the DP alt-mode link drops and renegotiates whenever ethernet and USB pull current.
- **Dock heat**: the link gets less stable as the temperature rises.
- **Cable and conversion path**: it is currently connected over HDMI, and **I have not yet
  tried a DP cable.**

> **Current status: unsolved.** The cause cannot be pinned down until a DP cable is tried.
> The cleanest way to separate the variables is to remove the dock and connect the MacBook to
> the monitor directly. If the symptom disappears it is the dock; if it persists it is the
> cable or the monitor.

---

## Wall 4. Waking a clamshell Mac remotely

This took the longest.

### Symptom

Cut the dock's power and the Mac sleeps; power it back on and no display appears. The keyboard
and mouse are attached to Windows, so there is no physical input device to wake the Mac with.

### Attempt 1: `acwake` (wake on power connection) — failed

```bash
pmset -g cap
```

The output had **no `acwake` entry**; this model does not support it.
It did have `womp` (Wake On Magic Packet).

### Attempt 2: wired Wake-on-LAN — failed

Magic packets sent from PowerShell got no response.

**Cause**: the Satechi dock's LAN port is a USB ethernet chip. When the Mac sleeps, USB power
drops and **the interface itself disappears.** On top of that, cutting the dock's power kills
the LAN port too, so there is no way for it to receive a packet.

### Attempt 3: Wi-Fi Wake-on-LAN — unreliable

It worked right after the Mac slept and stopped working later; results were erratic. A
sleeping Mac wakes briefly and periodically to run background work (dark wake), and a packet
that lands in that window succeeds while others fail. This is not something a setting can pin down.

### Attempt 4: SSH — worked

**An SSH connection wakes the Mac by itself.** No magic packet needed.

Enable remote login on the Mac:

```bash
sudo systemsetup -setremotelogin on
```

Then register a key from Windows, so the connection needs no password:

```powershell
ssh-keygen -t ed25519
Get-Content $env:USERPROFILE\.ssh\id_ed25519.pub | Set-Clipboard
```

Pass it to the Mac through Synergy's clipboard sharing, then in the Mac's terminal:

```bash
mkdir -p ~/.ssh && pbpaste >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys
```

> **Note**: existing files in the Windows `.ssh` folder do not interfere. Each direction needs
> a different file: the initiating side needs `id_ed25519` (the private key), the receiving
> side needs `authorized_keys` (the list of public keys).
> The opposite direction — Mac to Windows — is written up in
> [Eight Places I Got Stuck Setting Up SSH from a Mac to a Home Windows PC](/en/posts/mac-windows-ssh-setup/).
> On the Windows side there is a trap: for an administrator account the key goes in an
> entirely different file.

### Lighting the display took two commands

Even when SSH woke it, the Mac sometimes stayed in **dark wake** with no display. It turned
out that **two commands in succession** were needed to light it.

```powershell
ssh eddy@192.168.0.79 "caffeinate -u -t 5"
ssh eddy@192.168.0.79 "system_profiler SPDisplaysDataType"
```

The same command was useful for diagnosis.

```powershell
ssh eddy@192.168.0.79 "system_profiler SPDisplaysDataType | grep -E 'Online|Asleep'"
```

- `Display Asleep: Yes` → the Mac is awake and **only the display output** is asleep
- no `Online` line at all → the Mac does not even see the monitor

That distinction is what let me separate "network problem / Mac problem / monitor problem"
cleanly.

> **Current status: mostly working, with intermittent failures.**
> The clamshell display still occasionally fails to light. It may share a root cause with
> wall 3 (monitor flicker), so it needs rechecking after the DP cable swap.

---

## Wall 5. Switching monitor input and display mode together

### Goal

> The moment the dock powers on → Windows drops to a single monitor and the LG shows the Mac.
> Power the dock off → Windows returns to extended mode.

### The pieces

**Windows display mode**: `DisplaySwitch.exe`
(the same thing as the Win+P menu: `/extend`, `/external`, `/internal`, `/clone`)

**Monitor input switching**: ControlMyMonitor (DDC/CI, VCP code 60 = Input Select)

I found the LG UltraGear's input values by trial. In place of `<monitor serial>` below, put
your own monitor's serial number as reported by `ControlMyMonitor.exe /smonitors`.

```
ControlMyMonitor.exe /SetValue "<monitor serial>" 60 1   → Mac input
ControlMyMonitor.exe /SetValue "<monitor serial>" 60 3   → Windows input
```

I used the **serial number** as the identifier rather than an index like `\\.\DISPLAY2`,
because the index can change when the display mode changes.

> **Caution**: while the monitor is on the Mac input, entry 60 disappears from the listing
> entirely and cannot be read. **Writing still works**, though, so switching in both directions
> is possible. Other values such as brightness still read, so DDC communication itself is alive.

### Order mattered

**Once the LG is deactivated in Windows, DDC commands no longer reach it.** So the input switch
must always happen while the LG is still active in Windows.

- **To Mac mode**: switch input → wake the Mac → Windows single screen (mode change last)
- **To Windows mode**: sleep the Mac → Windows extend → (wait 3 s) → switch input

In AutoHotkey that becomes:

```autohotkey
; To Mac mode — input switch first, display mode change last
ToMacMode() {
    global dsw, lgMac, singleMode
    SetLgInput(lgMac)
    Sleep 2000
    WakeMac()
    Sleep 1000
    Run '"' dsw '" ' singleMode
}

; To Windows mode — revive the shared monitor with extend, then switch input
ToWinMode() {
    global dsw, lgWin
    SleepMac()
    Run '"' dsw '" /extend'
    Sleep 3000
    SetLgInput(lgWin)
}
```

### Detecting dock power

I tried sensing signal presence through the DDC `AE` value (Vertical Frequency), but it
returned 0 with the dock both on and off, so that failed.

Instead I judge it by **whether SSH connects.** With the dock off the Mac sleeps deeply and SSH
fails; with it on, the Mac answers. Check every 5 seconds with `ConnectTimeout=3` and act
**only at the moment the state changes.**

```autohotkey
macAlive := -1        ; -1 unknown, 1 on, 0 off
SetTimer CheckDock, 5000

CheckDock() {
    global ssh, mac, macAlive, autoConfirm

    code := RunWait('"' ssh '" -o BatchMode=yes -o ConnectTimeout=3 ' mac ' "true"',, "Hide")
    alive := (code = 0) ? 1 : 0

    if (alive = macAlive)
        return

    prev := macAlive
    macAlive := alive

    if (prev = -1)              ; the first check after startup only records
        return

    if (alive = 1)
        ToMacMode()
    else
        ToWinMode()
}
```

The full script is in `_code/synergy-mac-windows-setup/dual-control.ahk`
([view on GitHub](https://github.com/eddykim/eddykim.github.io/blob/main/_code/synergy-mac-windows-setup/dual-control.ahk)).
It also covers the Hangul/English pass-through, the manual switch shortcut (`Ctrl+Alt+D`), the
confirmation dialog before switching, and the `Sysnative` path handling that 32-bit AHK needs.
A shortcut in the `shell:startup` folder makes it run at boot.

> **Current status: half done.**
> Switching to Mac mode when the dock powers on works.
> The reverse — **automatically returning to Windows extended mode when the dock powers off —
> is still being refined.** For now I switch manually with `Ctrl+Alt+D`.

---

## Wall 6. The lock screen goes dark after one minute

### What I tried

`caffeinate -u -t 600` to hold it for ten minutes, which failed.
`pgrep -fl caffeinate` showed the process alive and well, and the screen went dark anyway.

### Cause

`caffeinate -u` sends activity signals to **the logged-in user session.** On the lock screen
the session is locked, so the signal never arrives. On top of that, the lock screen runs on
**its own short timer**, separate from the "turn display off" setting in System Settings.

I tried to raise it in settings, but the "Require password after screen saver begins or display
is turned off" item was **pinned to "immediately" and greyed out because of iPhone Mirroring's
automatic authentication.**

### Fix

Install **KeepingYouAwake.** Toggle it from the menu bar and that is the whole story. No
hidden `defaults` keys to poke, no security settings to weaken.

> **Lesson**: before digging into system settings, look for a well-made tool that already does it.

---

## Where it ended up

### Working

- One keyboard and mouse moving between the Mac and Windows
- One Hangul/English key switching input source on both operating systems
- `Ctrl+Alt+W` / `Ctrl+Alt+S` to wake and sleep the Mac
- `Ctrl+Alt+D` to switch manually between Mac mode and Windows extended mode
- Automatic switch to Mac mode when the dock powers on
- The LG serving as a Windows extended display while the Mac is idle
- KeepingYouAwake holding the lock screen on

### Still open

- **The clamshell display occasionally fails to light.** Cause unidentified.
- **Automatic reversion to Windows extended mode when the dock powers off needs work.**
  Manual for now.
- **Monitor flicker is not fully resolved.** To be rechecked after a DP cable swap.

Those three may well be entangled. The next step is to swap the cable and remove the hardware
variables first.

---

## What the detours taught me

**Read the error message precisely.**
"Serial key mismatch" was proof the network was already connected.
"Connection refused" meant the IP was right and no service was on the port.
When the message changes, the cause has changed.

**Separate symptoms by layer.**
"No display" covers three distinct states that look identical: ① the Mac is not waking,
② the Mac woke but its output is asleep, ③ the Mac is outputting but the monitor is not
receiving. One line — `system_profiler | grep -E 'Online|Asleep'` — split them apart.

**Check that the command exists on that OS.**
I used `setsid` at one point and lost a lot of time to it. It is a Linux command and macOS
does not have it. It failed silently, with no output, which made the cause slow to find.

**Prefer "change it into a structure that works" over "make the broken thing work."**
Instead of fighting to make Karabiner catch synthesized events, I moved the translation to
Windows. Instead of insisting on Wake-on-LAN for a sleeping Mac, I found a different path in
SSH. Looking for the side door beat pushing on the blocked one.

**Do not paper over a hardware problem with software.**
The flicker and the intermittent wake failures are most likely the cable or the dock, and I
spent a fair amount of time trying to solve them by adjusting script timing.

---

## References

- [Synergy purchase FAQ](https://symless.com/purchase/faq) — Symless; the source for personal tier prices
- [Differences between the business and personal editions](https://support.symless.com/hc/en-us/articles/33562900219665-What-are-the-differences-between-the-business-and-personal-editions) — Symless; computers per tier and whether TLS and hotkeys are included
- [Deskflow](https://github.com/deskflow/deskflow) — Synergy's upstream open-source project
- [Deskflow](https://wiki.archlinux.org/title/Deskflow) — Arch Wiki; the relationship with Synergy and the sponsorship
- [AutoHotkey v2 documentation](https://www.autohotkey.com/docs/v2/)
- [ControlMyMonitor](https://www.nirsoft.net/utils/control_my_monitor.html) — NirSoft; changes monitor input over DDC/CI
- [KeepingYouAwake](https://github.com/newmarcel/KeepingYouAwake) — holds a Mac awake
- [Karabiner-Elements](https://karabiner-elements.pqrs.org/) — Mac key remapping
