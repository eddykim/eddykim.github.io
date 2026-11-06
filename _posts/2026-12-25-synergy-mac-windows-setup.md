---
layout: post
title: "맥북과 윈도우 데스크탑을 키보드 하나로 쓰기까지 막힌 여섯 군데"
date: 2026-12-25 20:00:00 +0900
categories: [잡담]
page_id: synergy-mac-windows-setup
tags: [synergy, autohotkey, macos, windows, 클램쉘]
description: "맥북 클램쉘과 윈도우 데스크탑을 키보드·마우스 하나로 쓰고 모니터를 자동 전환하기까지 막힌 여섯 군데. 셋은 아직 미해결이다."
---

맥북(클램쉘, 사테치 큐브독)과 윈도우 데스크탑을 하나의 키보드·마우스로 제어하고,
모니터 두 대를 상황에 따라 자동으로 나눠 쓰는 환경을 만들었다.
단순해 보였지만 예상치 못한 벽이 여섯 개쯤 있었고, 그 과정을 정리한다.

> 이 글은 2026년 11월 기준 기록이다. 여섯 중 셋(모니터 깜박임, 클램쉘 디스플레이 점등,
> 자동 복귀)은 아직 해결되지 않았고, 각 절에 현재 상태를 적어 두었다.

## 환경

| 항목 | 내용 |
| --- | --- |
| 맥 | MacBook Pro (M4 Pro), 클램쉘 모드, 사테치 큐브독 연결 |
| 윈도우 | 데스크탑 (RTX 4080), **Synergy primary(서버)** |
| 모니터 1 | ASUS VA24DQLB — 윈도우 전용 |
| 모니터 2 | LG UltraGear — 윈도우·맥 공용, 입력 전환으로 사용 |
| 소프트웨어 | Synergy 3 Lite v3.6.3, AutoHotkey v2, ControlMyMonitor, KeepingYouAwake |

### 물리적 연결 구성

```
[윈도우 데스크탑] ──(DP/HDMI)──> [VA24DQLB]  … 윈도우 전용
        │
        └────────(HDMI)────────> [LG UltraGear] ─── 입력 3번
        │
        └────────(이더넷 직결)──> [사테치 큐브독] ── [MacBook Pro]
                                        │
                                        └──(HDMI)──> [LG UltraGear] ─── 입력 1번
```

핵심은 이렇다.

- **두 PC는 랜선으로 직결**되어 있다. 공유기를 거치지 않는 1:1 연결이고,
  맥 쪽은 사테치 독의 이더넷 포트를 쓴다.
- **데스크탑은 모니터 두 대를 모두 물고 있고, 맥북은 LG 한 대만 문다.**
  즉 LG 모니터는 두 컴퓨터에서 동시에 케이블이 꽂혀 있고, 입력 전환으로 주인이 바뀐다.
- 맥을 안 쓸 때는 LG를 윈도우의 확장 디스플레이로 활용하고,
  맥을 쓸 때는 LG를 맥 화면으로 넘기고 윈도우는 단일 모니터로 전환한다.

> 직결 랜선은 DHCP가 없어 양쪽 모두 `169.254.x.x`(APIPA) 주소를 받는다.
> Synergy 연결은 이 주소로 하고, 나중에 추가한 SSH 제어는 공유기를 거친 와이파이 주소를 쓴다.
> APIPA 주소는 재부팅 때 바뀔 수 있으므로 고정 IP로 바꿔두는 편이 안전하다.

---

## 문제 1. Synergy 연결 실패

### 증상

맥에서 윈도우 IP를 추가하려는데 **"Serial key mismatch"** 에러.

### 시도와 판단

- 네트워크 문제로 의심했으나, 이 에러가 떴다는 것 자체가 **상대 컴퓨터에 접속해 응답을 받아왔다는 증거**였다.
  연결이 안 됐다면 "host unreachable"이 떴을 것이다.
- 원인은 라이선스였다. 맥은 **Synergy 3 Max**, 윈도우는 다른 키를 쓰고 있었다.

### 해결

양쪽 모두 **Synergy 3 Lite**로 시리얼 키를 통일.
About → Update serial key → 같은 Symless 계정으로 로그인.

### 이어진 두 번째 실패

키를 맞추자 이번엔 **"Failed"**.
입력된 IP가 `192.254.214.194`였다. 첫 옥텟을 `169`로 쳐야 하는데 오타였다.

> **교훈**: 에러 메시지가 바뀌었다는 건 앞 단계가 해결됐다는 뜻이다.
> 같은 원인을 계속 의심하지 말고 새 메시지를 새로 읽자.

### 결과

`169.254.214.194`로 수정 후 정상 연결. 마우스·키보드 공유 성공.

---

## 곁가지: Synergy 라이선스 티어

시리얼 키를 맞추는 과정에서 든 의문이다.
어차피 양쪽 키가 같아야 하고 Synergy 자체 핫키 기능도 안 쓰는데, 티어를 올릴 이유가 있나?

Symless의 문서를 보면 개인용 티어는 Lite·Plus·Max가 각각 컴퓨터 3대·5대·15대를 지원하는
1회 결제 방식이고, **TLS/SSL 네트워크 암호화와 핫키는 Plus와 Max에만** 들어 있다.
가격은 2026년 11월 기준 Lite $15, Plus $19, Max $29다.
클립보드 공유는 모든 티어에 들어 있다.

즉 필자가 쓰는 **Lite는 암호화도 핫키도 없는 최소 티어**다.
컴퓨터 2대에 핫키를 안 쓴다면 티어를 올릴 이유는 없다는 결론이 나왔다.

무료 오픈소스 대안으로는 **Deskflow**가 있다. 흥미롭게도 Deskflow는 Synergy의
공식 업스트림 프로젝트이고, Synergy 측이 코드와 자금으로 후원하고 있다.
새로 시작하는 사람이라면 유료 결제 전에 이쪽을 먼저 시험해 볼 만하다.

---

## 곁가지: 클립보드 공유가 안 될 때

**증상**: 마우스·키보드 공유는 정상인데 윈도우에서 복사한 게 맥에서 붙여넣기가 안 됨.

**확인 순서**

1. Synergy → Settings → Copy & paste 토글 (기본 켜짐)
2. macOS 시스템 설정 → 개인정보 보호 및 보안 → **입력 모니터링**에 Synergy 추가

**원인**: 라이선스 변경으로 앱이 재설치되면서 권한이 초기화됐다.
키 입력 권한(손쉬운 사용)은 살아 있어 마우스·키보드는 멀쩡히 동작하는 탓에
클립보드만 안 되는 증상으로 나타나 원인을 짐작하기 어려웠다.

---

## 문제 2. 오른쪽 Alt 키가 맥에서 인식되지 않음

### 배경: 윈도우 주 사용자가 맥을 함께 쓸 때 겪는 문제

맥을 원래부터 쓰던 사람에게는 큰 이슈가 아니다.
맥 키보드에는 한영 전환 키가 따로 있고, Caps Lock으로도 전환되기 때문이다.

문제는 **필자처럼 윈도우를 주로 쓰다가 맥을 같이 쓰게 된 경우**다.
한글 키보드의 오른쪽 Alt는 윈도우에서 한/영 키(VK_HANGUL) 역할을 하는데,
맥에는 이 키에 대응하는 개념이 없다.

필자는 이미 맥 단독 환경에서 **Karabiner-Elements**로 우회하고 있었다.
macOS 한영 전환 단축키를 **F18**로 지정해두고,
Karabiner에서 **우측 Command → F18**로 매핑했다.
외부 키보드 기준으로 우측 Command 위치가 곧 오른쪽 Alt 자리이므로,
물리적으로는 "오른쪽 Alt로 한영 전환"이 되는 셈이었다.

### 증상

윈도우를 primary로 두고 Synergy로 연결하자, **이 Karabiner 매핑이 전혀 동작하지 않았다.**

### 왜 안 되는가

두 가지 이유가 겹쳐 있었다.

1. **Karabiner는 물리 HID 장치를 붙잡아(grab) 입력을 가로채는 방식**이다.
   Synergy가 맥에 넣는 키 입력은 물리 장치에서 온 게 아니라
   **소프트웨어가 합성해 주입하는 이벤트**라 Karabiner 레이어를 그냥 통과한다.
   즉 룰 자체가 발동하지 않는다.
   Karabiner-EventViewer를 띄워놓고 눌러보면 아무 이벤트도 안 찍힌다.

2. **Synergy의 수정자 키 매핑은 좌/우를 구분하지 못한다.**
   config에 `alt`라고 쓰면 좌우 Alt 전부에 적용되므로
   "오른쪽 Alt만 다르게"는 config 레벨에서 구조적으로 불가능하다.

### 해결: 변환을 윈도우 쪽으로 옮기기

macOS 한영 전환 단축키를 **F18에서 F16으로 변경**하고
(F18은 Karabiner가 이미 쓰고 있어 충돌을 피하려 바꿨다),
**윈도우에서 AutoHotkey로 F16을 직접 쏘는** 방식으로 전환했다.

Synergy는 이 F16을 그대로 맥에 전달하고,
macOS 시스템 단축키는 합성 이벤트에도 정상 반응한다. Karabiner를 거칠 필요가 없다.

```autohotkey
~vk15::Send "{Blind}{F16}"
```

- `vk15` = 한/영 키의 가상 키 코드
- `~` = 원래 키 기능도 통과시키므로 **윈도우의 한영 전환은 그대로 유지**되고 F16이 추가로 발생
- 커서가 윈도우에 있을 땐 F16이 아무 일도 안 하므로 부작용 없음

맥을 단독으로 쓸 때의 Karabiner 매핑(우측 Command → F18)은 그대로 남겨뒀다.
맥 자체 키보드로 쓸 때는 여전히 유효하기 때문이다.

### 결과

하나의 키로 양쪽 OS 한영 전환 성공.

---

## 문제 3. LG 모니터가 간헐적으로 깜박이거나 꺼짐

### 시도

Synergy의 `screenSaverSync = false`를 의심했으나 무관했다.
Synergy는 영상 출력에 전혀 관여하지 않는다.

### 부분적 원인

LG 모니터의 **초절전 전력 모드**였다.
대기 전력을 줄이려고 신호 감지 회로까지 꺼버려서,
독 전원이 다시 들어와도 모니터가 신호를 알아채지 못했다.
이 옵션을 끄자 "독을 켰는데 화면이 안 뜨는" 증상은 줄었다.

참고로 같은 메뉴의 **스마트 에너지 세이빙**은 화면 내용에 맞춰 밝기를 조절하는 기능이라
이 문제와 무관하다.

### 아직 남은 문제

사용 중 화면이 순간적으로 깜박이는 현상은 **여전히 남아 있다.**
의심되는 후보는 셋이다.

- **독의 전원 부족**: USB-C 독에 PD 충전기 없이 버스 파워만으로 쓰면
  이더넷·USB가 전류를 끌어갈 때 DP alt-mode 링크가 순간 끊겼다 재협상된다.
- **독의 발열**: 온도가 올라갈수록 링크가 불안정해진다.
- **케이블/변환 경로**: 현재 HDMI로 연결돼 있는데, **DP 케이블은 아직 시도해보지 않았다.**

> **현재 상태: 미해결.** DP 케이블로 바꿔보기 전까지는 원인을 특정할 수 없다.
> 가장 확실한 절분 방법은 독을 빼고 맥북에서 모니터로 직결해보는 것이다.
> 증상이 사라지면 독, 그대로면 케이블이나 모니터다.

---

## 문제 4. 클램쉘 맥을 원격으로 깨우기

가장 오래 걸린 부분이다.

### 증상

독 전원을 끄면 맥이 잠자기에 들어가고, 다시 켜도 화면이 안 뜬다.
키보드·마우스가 윈도우에 붙어 있으니 맥을 깨울 물리적 입력 장치가 없다.

### 시도 1: `acwake` (전원 연결 시 깨우기) — 실패

```bash
pmset -g cap
```

출력에 `acwake` 항목이 **없었다**. 이 모델은 미지원이다.
대신 `womp`(Wake On Magic Packet)는 있었다.

### 시도 2: 유선 Wake-on-LAN — 실패

PowerShell로 매직 패킷을 보냈지만 반응이 없었다.

**원인**: 사테치 독의 랜 포트는 USB 이더넷 칩이다.
맥이 잠들면 USB 전원이 내려가면서 **인터페이스 자체가 사라진다.**
게다가 독 전원을 끄면 랜 포트도 같이 죽으므로 패킷을 받을 방법이 없다.

### 시도 3: 와이파이 Wake-on-LAN — 불안정

잠든 직후에는 되고 시간이 지나면 안 되는 등 결과가 들쭉날쭉했다.
맥은 잠든 직후 백그라운드 작업을 위해 주기적으로 잠깐 깨어나는데(다크웨이크),
그 타이밍에 걸리면 성공하고 아니면 실패한다. 설정으로 고정할 수 있는 성질이 아니다.

### 시도 4: SSH — 성공

**SSH 접속 자체가 맥을 깨운다.** 매직 패킷이 필요 없었다.

맥에서 원격 로그인을 켜고:

```bash
sudo systemsetup -setremotelogin on
```

윈도우에서 키를 등록한다. 비밀번호 없이 접속하기 위한 단계다.

```powershell
ssh-keygen -t ed25519
Get-Content $env:USERPROFILE\.ssh\id_ed25519.pub | Set-Clipboard
```

Synergy 클립보드 공유로 맥에 넘긴 뒤, 맥 터미널에서:

```bash
mkdir -p ~/.ssh && pbpaste >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys
```

> **참고**: 윈도우 `.ssh` 폴더에 기존 파일이 있어도 무관하다.
> 방향마다 필요한 파일이 다르다. 거는 쪽은 `id_ed25519`(개인키),
> 받는 쪽은 `authorized_keys`(공개키 목록)다.
> 맥에서 윈도우로 거는 반대 방향 설정은
> [맥에서 집 윈도우 PC로 SSH 접속하기까지 막힌 여덟 군데](/posts/mac-windows-ssh-setup/)에 정리해 두었다.
> 윈도우 쪽은 관리자 계정이면 키를 넣는 파일이 아예 다르다는 함정이 있다.

### 디스플레이까지 켜려면 두 명령이 필요했다

SSH로 깨워도 **다크웨이크** 상태라 화면이 안 켜지는 경우가 있었다.
확인 결과 아래 **두 명령을 연달아** 실행해야 디스플레이가 켜졌다.

```powershell
ssh eddy@192.168.0.79 "caffeinate -u -t 5"
ssh eddy@192.168.0.79 "system_profiler SPDisplaysDataType"
```

상태 진단에도 이 명령이 유용했다.

```powershell
ssh eddy@192.168.0.79 "system_profiler SPDisplaysDataType | grep -E 'Online|Asleep'"
```

- `Display Asleep: Yes` → 맥은 깨어 있고 **디스플레이 출력만** 잠든 상태
- `Online` 줄 자체가 없음 → 맥이 모니터를 인식조차 못 하는 상태

이 구분 덕분에 "네트워크 문제 / 맥 문제 / 모니터 문제"를 정확히 나눌 수 있었다.

> **현재 상태: 대체로 동작하나 간헐적 실패 있음.**
> 클램쉘에서 디스플레이가 안 켜지는 경우가 아직 때때로 발생한다.
> 문제 3(모니터 깜박임)과 같은 뿌리일 가능성이 있어, DP 케이블 교체 후 재확인이 필요하다.

---

## 문제 5. 모니터 입력 + 디스플레이 모드 자동 전환

### 목표

> 독 전원을 켜는 순간 → 윈도우는 단일 모니터로 바뀌고, LG 모니터는 맥 화면을 띄운다.
> 독을 끄면 → 윈도우 확장 모드로 복귀.

### 구성 요소

**윈도우 디스플레이 모드**: `DisplaySwitch.exe`
(Win+P 메뉴와 동일. `/extend`, `/external`, `/internal`, `/clone`)

**모니터 입력 전환**: ControlMyMonitor (DDC/CI, VCP 코드 60 = Input Select)

LG UltraGear의 입력 값을 실측으로 찾았다. 아래 `<모니터 시리얼>` 자리에는
`ControlMyMonitor.exe /smonitors`로 조회한 자기 모니터의 시리얼 번호를 넣는다.

```
ControlMyMonitor.exe /SetValue "<모니터 시리얼>" 60 1   → 맥 입력
ControlMyMonitor.exe /SetValue "<모니터 시리얼>" 60 3   → 윈도우 입력
```

식별자로 `\\.\DISPLAY2` 같은 번호 대신 **시리얼 번호**를 썼다.
디스플레이 모드가 바뀌면 번호가 달라질 수 있기 때문이다.

> **주의**: 맥 입력 상태에서는 60번 항목이 목록에서 아예 사라진다(읽기 불가).
> 하지만 **쓰기는 정상 동작**해서 양방향 전환이 가능했다.
> 밝기 등 다른 값은 읽히므로 DDC 통신 자체는 살아 있는 상태였다.

### 순서가 중요했다

**LG가 윈도우에서 비활성화되면 DDC 명령이 닿지 않는다.**
그래서 입력 전환은 항상 LG가 윈도우에 활성화돼 있는 시점에 해야 한다.

- **맥 모드로**: 입력 전환 → 맥 깨우기 → 윈도우 단일 화면 (모드 변경이 마지막)
- **윈도우 모드로**: 맥 잠자기 → 윈도우 확장 → (3초 대기) → 입력 전환

AutoHotkey로 옮기면 이렇게 된다.

```autohotkey
; 맥 모드로 — 입력 전환을 먼저, 화면 모드 변경을 마지막에
ToMacMode() {
    global dsw, lgMac, singleMode
    SetLgInput(lgMac)
    Sleep 2000
    WakeMac()
    Sleep 1000
    Run '"' dsw '" ' singleMode
}

; 윈도우 모드로 — 확장으로 공용 모니터를 되살린 뒤 입력 전환
ToWinMode() {
    global dsw, lgWin
    SleepMac()
    Run '"' dsw '" /extend'
    Sleep 3000
    SetLgInput(lgWin)
}
```

### 독 전원 감지 방법

DDC의 `AE`(Vertical Frequency) 값으로 신호 유무를 감지하려 했으나
독 on/off 양쪽 모두 0이 나와 실패했다.

대신 **SSH 접속 성공 여부**로 판단했다.
독이 꺼져 있으면 맥이 깊이 잠들어 SSH가 실패하고, 켜지면 응답한다.
5초마다 `ConnectTimeout=3`으로 확인하고 **상태가 바뀌는 순간에만** 동작시킨다.

```autohotkey
macAlive := -1        ; -1 미확인, 1 켜짐, 0 꺼짐
SetTimer CheckDock, 5000

CheckDock() {
    global ssh, mac, macAlive, autoConfirm

    code := RunWait('"' ssh '" -o BatchMode=yes -o ConnectTimeout=3 ' mac ' "true"',, "Hide")
    alive := (code = 0) ? 1 : 0

    if (alive = macAlive)
        return

    prev := macAlive
    macAlive := alive

    if (prev = -1)              ; 시작 직후 첫 확인은 기록만
        return

    if (alive = 1)
        ToMacMode()
    else
        ToWinMode()
}
```

전체 스크립트는 `_code/synergy-mac-windows-setup/dual-control.ahk`에 있다
([GitHub에서 보기](https://github.com/eddykim/eddykim.github.io/blob/main/_code/synergy-mac-windows-setup/dual-control.ahk)).
한/영 전달, 수동 전환 단축키(`Ctrl+Alt+D`), 전환 전 확인창,
32비트 AHK의 `Sysnative` 경로 처리까지 들어 있다.
`shell:startup` 폴더에 바로가기를 넣어 부팅 시 자동 실행되게 했다.

> **현재 상태: 절반 완성.**
> 독을 켤 때 맥 모드로 넘어가는 건 동작한다.
> 반대로 **독을 껐을 때 윈도우 확장 모드로 자동 복귀하는 부분은 아직 다듬는 중**이다.
> 지금은 `Ctrl+Alt+D`로 수동 전환해서 쓰고 있다.

---

## 문제 6. 잠금 화면이 1분 만에 꺼짐

### 시도

`caffeinate -u -t 600`으로 10분 유지를 시도했으나 실패했다.
`pgrep -fl caffeinate`로 확인하면 프로세스는 멀쩡히 살아 있는데도 화면은 꺼졌다.

### 원인

`caffeinate -u`는 **로그인된 사용자 세션**에 활동 신호를 보내는 명령이다.
잠금 화면 상태에서는 세션이 잠겨 있어 신호가 닿지 않는다.
게다가 잠금 화면은 시스템 설정의 "디스플레이 끄기 시간"과 **별개의 짧은 타이머**로 동작한다.

설정에서 늘리려 했으나, "화면 보호기 시작 또는 디스플레이 꺼짐 후 암호 요구" 항목이
**iPhone 미러링 자동 인증 때문에 "즉시"로 고정**되어 회색 처리돼 있었다.

### 해결

**KeepingYouAwake** 설치. 메뉴바에서 켜고 끄는 것으로 끝이었다.
숨겨진 `defaults` 설정을 건드리거나 보안 설정을 낮출 필요가 없었다.

> **교훈**: 시스템 설정을 파고들기 전에 이미 잘 만들어진 도구가 있는지 먼저 찾아보자.

---

## 결과

### 완성된 것

- 하나의 키보드·마우스로 맥과 윈도우를 오가며 사용
- 한/영 키 하나로 양쪽 OS 한영 전환
- `Ctrl+Alt+W` / `Ctrl+Alt+S`로 맥 깨우기·재우기
- `Ctrl+Alt+D`로 맥 모드 ↔ 윈도우 확장 모드 수동 전환
- 독 전원을 켜면 자동으로 맥 모드 전환
- 맥을 안 쓸 때 LG 모니터를 윈도우 확장 디스플레이로 활용
- KeepingYouAwake로 잠금 화면 유지

### 남은 과제

- **클램쉘 디스플레이 점등이 때때로 실패한다.** 원인 미특정.
- **독을 껐을 때 윈도우 확장 모드로 자동 복귀하는 부분을 다듬어야 한다.** 현재는 수동 전환.
- **모니터 깜박임이 완전히 해결되지 않았다.** DP 케이블 교체 후 재검증 예정.

위 셋은 서로 얽혀 있을 가능성이 있다.
다음 단계는 케이블부터 바꿔서 하드웨어 변수를 먼저 제거하는 것이다.

---

## 삽질에서 얻은 것

**에러 메시지를 정확히 읽자.**
"Serial key mismatch"는 네트워크가 이미 연결됐다는 증거였다.
"Connection refused"는 IP는 맞고 포트에 서비스가 없다는 뜻이었다.
메시지가 바뀌면 원인도 바뀐 것이다.

**증상을 층별로 분리하자.**
"화면이 안 켜진다"는 ① 맥이 안 깨어남 ② 맥은 깨어났는데 출력이 잠듦
③ 맥은 출력하는데 모니터가 못 받음, 세 가지가 전부 같은 모습으로 나타난다.
`system_profiler | grep -E 'Online|Asleep'` 한 줄이 이 셋을 갈라줬다.

**명령어가 그 OS에 있는지 확인하자.**
중간에 `setsid`를 썼다가 한참 헤맸다. 리눅스 명령어고 macOS에는 없다.
아무 출력 없이 조용히 실패해서 원인을 찾는 데 오래 걸렸다.

**"안 되는 걸 되게 하기"보다 "되는 구조로 바꾸기".**
Karabiner로 합성 이벤트를 잡으려 애쓰는 대신 변환을 윈도우로 옮겼고,
잠든 맥을 Wake-on-LAN으로 깨우려 매달리는 대신 SSH 접속이라는 다른 경로를 찾았다.
막힌 방향을 계속 미는 것보다 옆문을 찾는 것이 빨랐다.

**소프트웨어로 하드웨어 문제를 덮으려 하지 말자.**
깜박임과 간헐적 점등 실패는 결국 케이블이나 독의 문제일 가능성이 높은데,
그걸 스크립트 타이밍 조정으로 해결하려다 시간을 꽤 썼다.

---

## 참고자료

- [Synergy 구매 FAQ](https://symless.com/purchase/faq) — Symless. 개인용 티어 가격의 출처
- [개인용·기업용 에디션 차이](https://support.symless.com/hc/en-us/articles/33562900219665-What-are-the-differences-between-the-business-and-personal-editions) — Symless. 티어별 컴퓨터 수와 TLS·핫키 포함 여부
- [Deskflow](https://github.com/deskflow/deskflow) — Synergy의 업스트림 오픈소스 프로젝트
- [Deskflow](https://wiki.archlinux.org/title/Deskflow) — Arch Wiki. Synergy와의 관계와 후원 내역
- [AutoHotkey v2 문서](https://www.autohotkey.com/docs/v2/)
- [ControlMyMonitor](https://www.nirsoft.net/utils/control_my_monitor.html) — NirSoft. DDC/CI로 모니터 입력을 바꾸는 도구
- [KeepingYouAwake](https://github.com/newmarcel/KeepingYouAwake) — 맥 잠금 화면 유지
- [Karabiner-Elements](https://karabiner-elements.pqrs.org/) — 맥 키 재매핑
