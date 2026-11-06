#Requires AutoHotkey v2.0
; ─────────────────────────────────────────────
; 맥 + 윈도우 듀얼 제어 스크립트
;
; 윈도우를 Synergy primary(서버)로 두고, 키보드·마우스 하나로 맥북(클램쉘)과
; 윈도우 데스크탑을 함께 쓰기 위한 보조 스크립트다. Synergy가 못 하는 세 가지,
; 즉 (1) 한/영 키를 맥에 전달하기, (2) 잠든 맥을 SSH로 깨우기,
; (3) 모니터 입력과 윈도우 디스플레이 모드를 함께 전환하기를 담당한다.
;
; 아래 네 값은 반드시 자기 환경 값으로 바꿔야 한다.
;   mac   : 맥의 사용자명@IP
;   cmm   : ControlMyMonitor.exe 경로
;   lg    : 공용 모니터의 시리얼 번호 (ControlMyMonitor /smonitors 로 확인)
;   lgMac / lgWin : 그 모니터의 입력 선택 값 (VCP 60의 실측값)
;
; 수정자 기호: ^ = Ctrl, ! = Alt, + = Shift, # = Win
; ~ = 원래 키 기능도 그대로 통과시킴
; vk15 = 한/영 키 (가상 키 코드 0x15)
; Run = 실행 후 바로 다음 줄 / RunWait = 끝날 때까지 대기
;
; 단축키
;   한/영      : 윈도우 한영 전환 + 맥 한영 전환(F16)
;   Ctrl+Alt+S : 맥 잠자기
;   Ctrl+Alt+W : 맥 깨우기
;   Ctrl+Alt+D : 모드 전환 (맥 모드 ↔ 윈도우 확장 모드)
;   자동       : 독 전원을 켜고 끄면 5초 내에 모드가 자동 전환
;
; 사전 준비
;   - 맥: 시스템 설정에서 원격 로그인(SSH)을 켜고, 윈도우의 공개키를
;     ~/.ssh/authorized_keys 에 등록해 비밀번호 없이 붙게 해 둔다.
;   - 맥: KeepingYouAwake 를 설치해 잠금 화면이 꺼지지 않게 한다.
;   - 윈도우: ControlMyMonitor(NirSoft)를 cmm 경로에 둔다.
;   - 이 파일의 바로가기를 shell:startup 에 넣으면 부팅 시 자동 실행된다.
; ─────────────────────────────────────────────

; ===== 준비 =====
; 32비트 AHK는 System32 접근이 다른 폴더로 돌려지므로 Sysnative를 먼저 시도
ssh := A_WinDir "\Sysnative\OpenSSH\ssh.exe"
if !FileExist(ssh)
    ssh := A_WinDir "\System32\OpenSSH\ssh.exe"

dsw := A_WinDir "\Sysnative\DisplaySwitch.exe"
if !FileExist(dsw)
    dsw := A_WinDir "\System32\DisplaySwitch.exe"

cmm := "C:\Tools\ControlMyMonitor\ControlMyMonitor.exe"
mac := "eddy@192.168.0.79"

lg    := "YOUR-MONITOR-SERIAL"   ; 공용 모니터 시리얼 번호
lgMac := 1                       ; 맥(독) 입력
lgWin := 3                       ; 윈도우 입력

; 맥 모드일 때 윈도우가 남길 화면
; /external = "두 번째 화면만", /internal = "PC 화면만"
singleMode := "/external"

; ===== 공용 함수 =====

; 맥 깨우기 — 두 명령을 연달아 실행해야 디스플레이가 켜짐 (실측 확인)
;   화면 유지는 맥의 KeepingYouAwake가 담당
WakeMac() {
    global ssh, mac
    RunWait '"' ssh '" -o BatchMode=yes -o ConnectTimeout=15 ' mac ' "caffeinate -u -t 5"',, "Hide"
    RunWait '"' ssh '" -o BatchMode=yes -o ConnectTimeout=15 ' mac ' "system_profiler SPDisplaysDataType"',, "Hide"
}

SleepMac() {
    global ssh, mac
    Run '"' ssh '" -o BatchMode=yes ' mac ' "pmset sleepnow"',, "Hide"
}

; 공용 모니터 입력 전환 (VCP 60 = Input Select)
SetLgInput(value) {
    global cmm, lg
    RunWait '"' cmm '" /SetValue "' lg '" 60 ' value,, "Hide"
}

; 맥 모드로 — 입력 전환을 먼저, 화면 모드 변경을 마지막에
; (공용 모니터가 윈도우에서 비활성화되면 DDC 명령이 닿지 않음)
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

; ===== 단축키 =====

~vk15::Send "{Blind}{F16}"      ; 한/영: 윈도우 전환 유지 + 맥으로 F16
^!s::SleepMac()
^!w::WakeMac()

^!d:: {
    global macAlive
    if MonitorGetCount() >= 2 {
        ToMacMode()
        macAlive := 1           ; 자동 감지와 상태를 맞춰둠
    } else {
        ToWinMode()
        macAlive := 0
    }
}

; ===== 독 전원 자동 감지 =====
; SSH 접속 성공 = 독 켜짐, 실패 = 독 꺼짐
; 상태가 바뀌는 순간에만 동작

macAlive := -1        ; -1 미확인, 1 켜짐, 0 꺼짐
autoConfirm := true   ; 전환 전 확인창 (익숙해지면 false)

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

    if (alive = 1) {
        if (autoConfirm && MsgBox("맥 모드로 전환할까요?",, "YesNo") != "Yes")
            return
        ToMacMode()
    } else {
        if (autoConfirm && MsgBox("윈도우 확장 모드로 전환할까요?",, "YesNo") != "Yes")
            return
        ToWinMode()
    }
}
