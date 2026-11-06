---
layout: post
title: "맥에서 집 윈도우 PC로 SSH 접속하기까지 막힌 여덟 군데"
date: 2026-12-23 20:00:00 +0900
categories: [잡담]
page_id: mac-windows-ssh-setup
tags: [ssh, windows, openssh, ddns, 포트포워딩]
description: "맥에서 집 윈도우 PC로 SSH 접속 환경을 만들며 걸린 여덟 가지 함정과, 그때마다 엉뚱한 곳을 보고 있었다는 기록."
---

집에 있는 윈도우 데스크탑에 맥에서 SSH로 접속해 작업하는 환경을 만들었다.
"OpenSSH 설치하고 포트 열면 끝"일 줄 알았는데, 실제로는 예상 못 한 지점에서 계속 막혔다.
이 글은 성공한 설정 방법보다 **어디서 어떻게 막혔고 왜 그랬는지**에 무게를 둔 기록이다.

윈도우에 OpenSSH 서버를 올리는 기본 과정 자체는
[이 글](https://blog.naver.com/twilight_teatime/223668641470)과
[마이크로소프트 설치 문서](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_install_firstuse)를
참고했다. 여기서는 그 이후의 이야기를 다룬다.

> 이 글의 IP·도메인·포트는 전부 예시값이다. 본인 환경 값으로 바꿔서 읽어야 한다.

---

## 최종 구성

```
[맥]  ──인터넷──▶  공유기(공인 IP)  ──포트포워딩──▶  [윈도우 PC]
                   외부 포트 1927        192.168.0.59:4131
                         │
                  DDNS(No-IP)로
                myhost.ddns.net 유지
```

- 윈도우: OpenSSH Server, 포트 `4131` (기본 22에서 변경)
- 공유기: 외부 `1927` → 내부 `192.168.0.59:4131` 포트포워딩
- DDNS: 유동 공인 IP를 `myhost.ddns.net`으로 고정
- 인증: ed25519 공개키
- 클라이언트: 맥 터미널 + VS Code Remote-SSH

---

## 걸린 함정 여덟 개

먼저 결론부터. 필자가 걸린 함정은 총 여덟 개였다.

| # | 증상 | 진짜 원인 |
|---|---|---|
| 1 | 공개키를 넣어도 비밀번호를 계속 물음 | 관리자 계정은 `authorized_keys`를 안 읽는다 |
| 2 | 키를 등록했는데 인증이 안 됨 | 안내문의 **자리표시자를 그대로** 붙여넣음 |
| 3 | `pbcopy: command not found` | SSH 세션 안에서 맥 명령을 실행 |
| 4 | 예전 설정이 전부 무효 | 공유기가 바뀌어 IP 대역이 통째로 변경 |
| 5 | DynDNS를 고를까 고민 | 무료 플랜이 이미 폐지되어 선택지가 아님 |
| 6 | 공유기 DDNS 인증 실패 | DDNS Key 방식은 도메인 칸에 다른 값을 넣어야 함 |
| 7 | 포워딩했는데 timeout | 집 안에서 테스트해서 NAT 루프백에 걸림 |
| 8 | 터미널은 되는데 VS Code만 실패 | 입력창에 `ssh` 명령어를 통째로 넣음 |

---

## 함정 1. 관리자 계정은 `authorized_keys`를 읽지 않는다

가장 먼저, 그리고 가장 오래 헤맨 지점이다.

공개키를 `C:\Users\사용자\.ssh\authorized_keys`에 넣었는데 계속 비밀번호를 물었다.
권한도 맞췄고 키 형식도 멀쩡한데 왜 안 되나 싶었는데, `sshd_config` 맨 아래에 답이 있었다.

```
Match Group administrators
       AuthorizedKeysFile __PROGRAMDATA__/ssh/administrators_authorized_keys
```

윈도우용 OpenSSH는 **Administrators 그룹에 속한 계정에만 별도의 공용 파일**을 본다.
필자의 계정이 관리자였으므로, 홈 디렉토리의 `authorized_keys`는 처음부터 쳐다보지도 않고 있었던 것이다.

등록해야 할 진짜 위치는 여기다.

```
C:\ProgramData\ssh\administrators_authorized_keys
```

자기 계정이 관리자인지는 이렇게 확인한다.

```powershell
Get-LocalGroupMember -Group Administrators
```

권한도 중요하다. 이 파일은 **Administrators와 SYSTEM만** 접근 가능해야 하고,
아니면 sshd가 조용히 무시한다. 에러도 안 뱉는다.

```powershell
icacls C:\ProgramData\ssh\administrators_authorized_keys /inheritance:r /grant "Administrators:F" /grant "SYSTEM:F"
Restart-Service sshd
```

이 동작은 마이크로소프트의 [OpenSSH 키 관리 문서](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_keymanagement)에
명시되어 있다. 문서를 먼저 읽었다면 아낄 수 있던 시간이다.

> **교훈:** 윈도우에서 공개키 인증이 안 되면 파일 위치부터 의심하자.
> 리눅스 경험만 있으면 절대 예상 못 하는 동작이다.

---

## 함정 2. 자리표시자를 그대로 붙여넣었다

안내에 이런 명령어가 있었다.

```powershell
Add-Content -Path C:\ProgramData\ssh\administrators_authorized_keys -Value 'ssh-ed25519 AAAA...복사한내용... eddy@mac' -Encoding ascii
```

`AAAA...복사한내용...` 자리에 진짜 키를 넣어야 하는데, **통째로 그대로 실행**해버렸다.
당연히 인증은 안 됐고, 파일에는 쓰레기 한 줄이 남았다.

재밌는 건 이걸 **파일 크기 변화만으로** 잡아냈다는 점이다.

```
등록 전: 749 바이트
등록 후: 787 바이트   → +38
```

ed25519 공개키 한 줄은 보통 90바이트 안팎이다. 38바이트는 너무 작다.
그리고 `-Encoding ascii`로 저장하면 한글이 `?` 한 글자씩으로 치환되는데,
`'ssh-ed25519 AAAA...복사한내용... eddy@mac'`를 그렇게 계산하면 정확히 38바이트가 나온다.
`ssh-ed25519`와 공백이 12, `AAAA...복사한내용...`이 15(한글 다섯 자가 `?` 다섯 개),
공백 1, `eddy@mac`이 8, 줄바꿈 `CRLF`가 2다.

정리는 이렇게 했다.

```powershell
$p = 'C:\ProgramData\ssh\administrators_authorized_keys'
(Get-Content $p) | Where-Object { $_ -notmatch '\.\.\.' } | Set-Content $p -Encoding ascii
```

> **교훈:** 공개키를 등록했으면 **파일 크기부터 확인하자.**
> ed25519면 +90바이트 안팎, RSA 4096이면 +750바이트 안팎이 정상이다.
> sshd는 잘못된 키 줄을 만나도 조용히 건너뛰기 때문에 에러로는 알 수가 없다.

---

## 함정 3. SSH 세션 안에서 맥 명령을 실행했다

맥에서 키를 만드는 순서는 이랬다.

```bash
ssh-keygen -t ed25519 -C "eddy@mac"
pbcopy < ~/.ssh/id_ed25519.pub
ssh -p 4131 eddy@192.168.0.59
```

그런데 필자는 **SSH 접속을 먼저 하고** 저 명령들을 실행했다. 결과는 이렇다.

- `ssh-keygen` → 윈도우에 키가 생성됨 (맥이 아니라!)
- `pbcopy` → macOS 전용 명령이라 윈도우에서 "없는 명령" 에러

SSH로 접속한 뒤에는 터미널 창은 맥 것이지만 **키보드 입력은 전부 원격 PC로 간다.**
머리로는 알면서도 창 모양이 똑같아서 계속 헷갈렸다.

지금은 헷갈릴 때마다 이걸 친다.

```bash
uname -s
```

- `Darwin` → 맥
- 에러 또는 다른 출력 → 원격(윈도우)

**터미널 탭을 두 개 열어 하나는 로컬 전용, 하나는 SSH 전용으로 고정**하는 것이 근본 해결책이었다.

---

## 함정 4. 공유기를 바꾸면 예전 메모는 전부 쓰레기가 된다

예전에 적어둔 메모에는 이렇게 있었다.

```
203.0.113.17:1927 -> 192.168.55.27:4131
```

그런데 실제로 확인해보니 현재 값은 전부 달랐다.

```
공인 IP     203.0.113.17    →  203.0.113.88
내부 IP     192.168.55.27   →  192.168.0.59
게이트웨이  192.168.55.1    →  192.168.0.1
```

이사나 공유기 교체로 **대역이 통째로 바뀌면** 포트포워딩 규칙은 새 공유기에 없다.
설정을 시작하기 전에 현재 값부터 다시 확인해야 헛수고를 안 한다.

```powershell
Get-NetIPAddress -AddressFamily IPv4 | Where-Object { $_.PrefixOrigin -eq 'Dhcp' }
Get-NetRoute -DestinationPrefix '0.0.0.0/0' | Select-Object NextHop
(Invoke-RestMethod 'https://api.ipify.org?format=json').ip
```

### 반드시 먼저 할 것: CGNAT 검사

DDNS와 포트포워딩 설정에 시간을 쏟기 전에 **이것부터** 확인해야 한다.

공유기 관리자 페이지의 WAN IP와, 위에서 조회한 공인 IP를 비교한다.

| 공유기 WAN IP | 판정 |
|---|---|
| 공인 IP와 동일 | 정상. 포트포워딩 가능 |
| `100.64.x.x` ~ `100.127.x.x` | CGNAT. 포트포워딩 불가 |
| `10.x.x.x`, `192.168.x.x`, `172.16~31.x.x` | 이중 NAT. 앞단에 다른 공유기 존재 |

아래 두 경우면 포트포워딩도 DDNS도 작동하지 않는다.
통신사에 공인 IP를 요청하거나, **Tailscale 같은 오버레이 VPN으로 전환**하는 것이 답이다.
솔직히 공유기를 건드릴 필요도 없고 포트를 인터넷에 노출하지도 않으니 그쪽이 더 편하다.
포트포워딩을 고집할 이유는 "추가 소프트웨어를 안 깔고 싶다" 정도뿐이다.

필자의 경우는 다행히 WAN IP가 공인 IP와 일치해서 그대로 진행했다.

---

## 함정 5. DynDNS는 이미 선택지가 아니다

공유기 DDNS 설정에 제공업체가 `NO-IP`와 `DynDNS` 둘뿐이었다.

그런데 **DynDNS의 무료 호스트네임 서비스는 2014년 5월에 끝났다.**
Dyn의 공지에 따르면 무료 계정은 그 시점에 비활성화되고 유료 Remote Access(당시 DynDNS Pro)로
넘어가야 했다. 지금은 Oracle이 그 서비스를 운영하는 유료 상품이다.
공유기 펌웨어에 항목이 남아 있는 건 그냥 오래된 목록일 뿐이다.

그래서 **No-IP**를 선택했다. 무료 플랜으로 호스트 3개까지 쓸 수 있다.

### No-IP 호스트 생성

`My No-IP → Dynamic DNS → Hostnames → Create Hostname`

| 항목 | 값 |
|---|---|
| Type | `A` |
| Host | 원하는 이름 (전 세계에서 유일해야 함) |
| Domain | `ddns.net` 등 무료 도메인 |
| IPv4 | 자동으로 채워진 현재 공인 IP |
| TTL | Auto |
| Wildcard | 해제 (유료 기능) |
| **Enable Dynamic DNS** | **체크** |

마지막 항목이 핵심이다. 체크하지 않으면 **고정 A 레코드**가 만들어져서,
IP가 바뀌어도 공유기가 갱신하지 못한다. DDNS를 쓰는 이유 자체가 사라진다.

### 무료 플랜의 함정

무료 No-IP는 **30일마다 확인 메일의 링크를 눌러야** 한다. 안 누르면 호스트가 삭제된다.
잘 쓰던 게 어느 날 갑자기 안 되면 십중팔구 이것이다.

### 도메인 선택

`ddns.net`, `hopto.org`처럼 널리 쓰이는 무료 DDNS 도메인은
기업 방화벽이나 공공 와이파이 DNS 필터에 차단된 사례가 보고된다.
필자가 여러 망에서 측정해 본 것은 아니니 일반화할 수는 없지만,
회사에서 접속할 일이 많다면 접속이 막힐 수 있다는 가능성은 염두에 둘 만하다.
No-IP는 무료 플랜에서도 여러 도메인을 고를 수 있으니, 막히면 다른 도메인으로 호스트를 다시 만들면 된다.

---

## 함정 6. DDNS Key 방식은 도메인 칸에 다른 값을 넣는다

`Enable Dynamic DNS`를 체크하면 **DDNS Key**를 발급받는 화면이 나온다.
공유기에 No-IP 계정 비밀번호를 저장하지 않아도 되는 방식이라 이쪽이 더 안전하다.

여기서 헷갈린 지점. 발급 화면에 이런 안내가 있다.

> Use Hostname for Integrated DDNS Devices: **all.ddnskey.com**

즉 공유기에 입력할 값은 이렇게 된다.

| 공유기 필드 | 입력값 |
|---|---|
| 서비스 제공 업체 | `NO-IP` |
| 사용자이름 | DDNS Key 사용자명 (`xxxxxxx@ddnskey.com` 형식) |
| 비밀번호 | DDNS Key 비밀번호 |
| 도메인 이름 | **`all.ddnskey.com`** ← 자기 호스트네임이 아니다 |
| WAN IP 바인딩 | **해제** |

`all.ddnskey.com`은 "이 키에 연결된 모든 호스트네임을 갱신하라"는 뜻의 특수 주소다.
공유기가 여기로 신호를 보내면 No-IP가 알아서 `myhost.ddns.net`을 갱신한다.
No-IP의 [공유기 DDNS 설정 문서](https://www.noip.com/support/knowledgebase/how-to-configure-ddns-in-router)도
도메인 칸에 이 값을 넣으라고 안내한다.

**접속할 때 쓰는 주소는 여전히 `myhost.ddns.net`이다.**
`all.ddnskey.com`은 공유기 설정 칸에만 들어가는 값이다. 이 구분이 제일 헷갈렸다.

> 실제로는 호스트네임(`myhost.ddns.net`)을 직접 넣어도 동작했다.
> DDNS Key가 그 호스트네임에 할당돼 있으면 둘 다 된다.
> 인증 실패가 뜨면 `all.ddnskey.com`으로 바꿔보면 된다.

**비밀번호는 발급 시 딱 한 번만 표시된다.** 메모장을 먼저 열어두고 생성 버튼을 누르자.
재발급하면 공유기 설정을 처음부터 다시 해야 한다.

### WAN IP 바인딩은 켜지 말 것

이 옵션을 켜면 공유기가 **자기 WAN 인터페이스 IP를 그대로** No-IP에 보고한다.
CGNAT 환경이라면 사설 IP가 DNS에 등록되면서 **멀쩡하던 레코드까지 망가진다.**

해제해두면 No-IP가 접속해온 패킷의 출발지 IP, 즉 진짜 공인 IP를 자동 감지한다.
일반 가정용 공유기에서는 해제가 정답이다.

---

## 함정 7. 집 안에서 테스트하면 실패한다

포트포워딩을 설정했다. 아래는 TP-Link 공유기의 메뉴 경로다.

```
고급 → NAT 포워딩 → 가상 서버 → 추가
  외부 포트  1927
  내부 IP    192.168.0.59
  내부 포트  4131
  프로토콜   TCP
```

> 같은 메뉴의 **포트 트리거링**과 혼동하기 쉽다. **가상 서버(Virtual Servers)** 가 맞다.

그런데 집 와이파이에 연결된 맥에서 테스트하니 timeout이 났다.
설정이 틀린 줄 알고 한참 규칙을 뒤졌는데, 원인은 **NAT 루프백(헤어핀 NAT)** 이었다.

내부망에서 자기 공인 IP로 나갔다가 다시 들어오는 경로인데, **지원하지 않는 공유기가 많다.**
포트포워딩이 완벽해도 집 안에서는 실패할 수 있다.

**반드시 휴대폰 핫스팟 등 외부망에 연결해서 테스트해야 한다.** 이때 맥의 와이파이는 끈다.
포트가 실제로 열렸는지만 확인하려면 [canyouseeme.org](https://canyouseeme.org/) 같은
외부 포트 점검 서비스로도 된다.

그리고 내부 IP 고정도 같이 해두자. DHCP로 IP가 바뀌면 포워딩이 조용히 깨진다.

```
고급 → 네트워크 → DHCP 서버 → 주소 예약
```

---

## 함정 8. 터미널은 되는데 VS Code만 안 된다

터미널에서는 접속이 됐는데 VS Code Remote-SSH만 이런 에러를 뱉었다.

```
WARN: $PLATFORM is undefined in installation script output. Errors may be dropped.
Failed to parse remote port from server output
Exec server for ssh-remote+7b22686f73744e616d65223a226d79686f73742e64646e732e6e6574222c...
```

에러 메시지가 불친절해 보이지만, **뒤의 16진수는 그냥 JSON을 hex로 인코딩한 것**이다. 디코딩하면 이렇다.

```json
{"hostName":"myhost.ddns.net","user":"ssh -p 1927 eddy"}
```

범인이 바로 나왔다. 입력창에 `ssh -p 1927 eddy@myhost.ddns.net`을 통째로 붙여넣었더니
**`ssh -p 1927 eddy` 전체가 사용자 이름으로** 들어갔다. 포트는 기본값 22로 시도했으니 될 리가 없다.

디코딩은 이렇게 한다.

```bash
echo '7b22686f73744e616d65...' | xxd -r -p
```

> **교훈:** Remote-SSH 에러의 16진수는 반드시 디코딩해보자. 원인이 그 안에 다 적혀 있다.

### 해결: `~/.ssh/config`에 적어둔다

입력창에 매번 치는 방식은 이런 실수가 반복된다. **맥의** `~/.ssh/config`에 한 번 적어두면 끝난다.

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

이제 터미널에서도 짧아지고,

```bash
ssh winpc        # 집에서
ssh winpc-ext    # 밖에서
```

VS Code에서는 `Cmd+Shift+P → Remote-SSH: Connect to Host...` 를 열면
**`winpc`, `winpc-ext`가 목록에 그냥 뜬다.** 더 이상 아무것도 입력할 필요가 없다.

> 설정 파일은 VS Code에서 `Cmd+Shift+P → Remote-SSH: Open SSH Configuration File...` 로 바로 열 수 있다.
> 첫 접속은 1~2분 걸린다. VS Code가 원격 PC에 서버 컴포넌트를 설치하기 때문이다.

---

## 마무리 체크리스트

순서대로 하면 위의 삽질은 대부분 피할 수 있다.

1. **CGNAT 검사** — 공유기 WAN IP와 공인 IP 비교. 다르면 여기서 멈추고 Tailscale로 전환
2. 윈도우에 OpenSSH Server 설치, 포트 변경 시 **방화벽 규칙도 그 포트로** 추가
3. 관리자 계정이면 공개키는 `C:\ProgramData\ssh\administrators_authorized_keys`
4. 키 등록 후 **파일 크기 증가분 확인** (ed25519 기준 +90바이트 안팎)
5. **같은 공유기 안에서** 먼저 접속 성공시키기
6. No-IP 호스트 생성(`Enable Dynamic DNS` 체크) 후 공유기 DDNS 설정(WAN IP 바인딩 해제)
7. 포트포워딩 + DHCP 주소 예약
8. **핫스팟으로 전환 후** 외부 접속 테스트
9. 맥 `~/.ssh/config` 작성
10. VS Code Remote-SSH는 목록에서 선택

## 진단 명령 모음

윈도우 쪽 상태 확인:

```powershell
Get-Service sshd                                             # 서비스 상태
Get-NetTCPConnection -State Listen -LocalPort 4131           # 리스닝 확인
Get-Item C:\ProgramData\ssh\administrators_authorized_keys   # 키 파일 크기
Resolve-DnsName myhost.ddns.net -Type A -DnsOnly             # DDNS 반영 확인
(Invoke-RestMethod 'https://api.ipify.org?format=json').ip   # 현재 공인 IP
```

맥 쪽:

```bash
ssh -v -p 1927 eddy@myhost.ddns.net    # 어디서 멈추는지 추적
dig +short myhost.ddns.net             # DNS 조회
uname -s                               # 로컬인지 원격인지
```

`ssh -v` 결과로 보는 원인:

| 메시지 | 원인 |
|---|---|
| `Operation timed out` | 포트포워딩 누락, 또는 NAT 루프백 |
| `Connection refused` | 포워딩은 됐으나 내부 포트가 틀림 |
| `Could not resolve hostname` | DDNS 미반영 또는 오타 |
| `Permission denied (publickey,password)` | **네트워크는 성공.** 키 등록 문제만 남음 |

마지막 줄이 뜨면 절반은 성공한 것이다. 네트워크 경로는 뚫렸고 인증만 남았다는 뜻이다.

---

## 보안 정리

접속이 되고 나면 마무리로 해둘 것들.

**1. 비밀번호 인증 끄기**

공개키 접속이 확실히 되는 걸 확인한 **뒤에**, 접속된 창 하나는 열어둔 채로 바꾼다.
`sshd_config`는 고치기 전에 사본을 떠 두는 편이 안전하다.

```
# C:\ProgramData\ssh\sshd_config
PasswordAuthentication no
```

포트를 인터넷에 여는 순간 전 세계에서 스캔이 들어온다. 비밀번호 인증은 닫아두는 게 좋다.

**2. 자격증명을 평문으로 남기지 않기**

설정하면서 만든 메모 파일에 키 passphrase나 DDNS Key 비밀번호를 적어두기 쉽다.
PC가 뚫리면 같이 털리니 비밀번호 관리자로 옮기자.

**3. 기본 셸 바꾸기 (선택)**

윈도우 OpenSSH의 기본 셸은 `cmd.exe`다. PowerShell이 편하면 이렇게 바꾼다.

```powershell
New-ItemProperty -Path "HKLM:\SOFTWARE\OpenSSH" -Name DefaultShell `
  -Value "C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe" -PropertyType String -Force
```

단, VS Code Remote-SSH는 `cmd.exe`에서 가장 안정적으로 동작한다.
Remote-SSH를 주로 쓴다면 그냥 두는 편이 낫다.

---

## 돌아보며

여덟 함정 중 **여섯 개가 "설정값이 틀렸다"가 아니라 "엉뚱한 곳을 보고 있었다"** 였다.

- 키를 엉뚱한 파일에 넣고 (함정 1)
- 엉뚱한 내용을 넣고 (함정 2)
- 엉뚱한 기계에서 명령을 치고 (함정 3)
- 엉뚱한 네트워크의 옛날 값을 보고 (함정 4)
- 엉뚱한 칸에 호스트네임을 넣고 (함정 6)
- 엉뚱한 위치에서 테스트했다 (함정 7)

SSH 설정이 어려운 게 아니라, **지금 어느 기계의 어느 파일을 보고 있는지**를
놓치기 쉬운 구조라는 게 진짜 난이도였다.

그래서 가장 도움이 된 습관 두 가지로 마무리한다.

- **의심스러우면 `uname -s`** — 지금 로컬인가 원격인가
- **바꿨으면 크기와 조회 결과를 확인** — 에러 없이 조용히 무시되는 설정이 많다

---

## 참고자료

- [Windows에서 OpenSSH 설치하기](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_install_firstuse) — Microsoft Learn
- [OpenSSH 키 관리](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_keymanagement) — Microsoft Learn. `administrators_authorized_keys` 동작의 출처
- [윈도우 SSH 서버 구축](https://blog.naver.com/twilight_teatime/223668641470) — 기본 설치 과정 참고
- [OpenSSH 서버 포트 번호 변경](https://www.lainyzine.com/ko/article/how-to-change-port-number-of-openssh-server-on-windosw-10/) — lainyzine
- [DynDNS 무료 사용자 종료 FAQ](https://help.dyn.com/free-user-discontinuation-faqs/) — Dyn. 무료 플랜 종료 시점의 출처
- [공유기에 DDNS 설정하기](https://www.noip.com/support/knowledgebase/how-to-configure-ddns-in-router) — No-IP
- [DDNS Key 만들고 쓰기](https://www.noip.com/support/knowledgebase/how-to-setup-and-use-a-ddns-key) — No-IP
- [TP-Link 공유기 DDNS 설정](https://www.noip.com/support/knowledgebase/setting-ddns-tp-link-router) — No-IP
- [Remote Development using SSH](https://code.visualstudio.com/docs/remote/ssh) — VS Code 문서
- [Open Port Check Tool](https://canyouseeme.org/) — 외부에서 포트가 열렸는지 확인
