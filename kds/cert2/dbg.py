import sys,runpy
sys.argv=['surf_run.py','floorc2','0.9','0','0','0','1','0.1','1','1']
import builtins
src=open('surf_run.py').read().replace("r,deep=bb(BB,rule,4,maxdep=44)","b=BB\nprint({k:(mn(b[k])>=0, mx(b[k])<=0) for k in b})\nr,deep=bb(BB,rule,4,maxdep=44)")
exec(compile(src,'x','exec'))
