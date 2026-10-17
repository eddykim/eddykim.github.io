"""그림·표 번호 뒤 조사가 숫자 읽기와 맞는지 본다 (그림3은 → 그림3은, 그림2은 → 그림2는).
figures.py renumber는 번호만 바꾸고 조사는 고치지 않으므로, 번호를 바꾼 뒤 돌린다.
사용: python3 josa_check.py "<한국어판 경로 또는 glob>"   (-en.md는 건너뛴다)"""
import re,sys,glob
# 숫자 읽기의 받침: 0영 1일 2이 3삼 4사 5오 6육 7칠 8팔 9구, 10의 배수는 '십'(ㅂ)
def batchim(num):
    n=int(num)
    if n%10==0 and n>0: return 'ㅂ'
    d=n%10
    return {0:'ㅇ',1:'ㄹ',2:None,3:'ㅁ',4:None,5:None,6:'ㄱ',7:'ㄹ',8:'ㄹ',9:None}[d]
pairs={'은':('은','는'),'는':('은','는'),'이':('이','가'),'가':('이','가'),'을':('을','를'),'를':('을','를'),'과':('과','와'),'와':('과','와')}
pat=re.compile(r'(그림|표)\s?(\d+)(\([a-z]\))?(은|는|이|가|을|를|과|와|으로|로)(?![가-힣])')
bad=0
for f in sorted(glob.glob(sys.argv[1])):
    if f.endswith('-en.md'): continue
    for i,line in enumerate(open(f,encoding='utf-8'),1):
        for m in pat.finditer(line):
            b=None if m.group(3) else batchim(m.group(2)); j=m.group(4)  # 패널 기호 a~d는 받침 없음
            if j in ('으로','로'):
                want='로' if (b is None or b=='ㄹ') else '으로'
            else:
                want=pairs[j][0] if b else pairs[j][1]
            if j!=want:
                bad+=1; print(f"{f}:{i}: {m.group(0)} → {m.group(1)}{m.group(2)}{m.group(3) or ''}{want}")
print("총",bad,"건")
