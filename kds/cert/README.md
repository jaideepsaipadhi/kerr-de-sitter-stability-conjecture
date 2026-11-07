# Event-lemma certification (computer-assisted proof)
- verifier.py: rigorous point verifier (interval arithmetic, bisection in z).
- boxes2.py: box verifier.  Center-fixed multiplier; (A,B) via corners; adaptive parameter splitting with a box cap.
- driver.py: parallel tiling with a resumable log (cert_log.txt).
    python3 driver.py small     # quick
    python3 driver.py core      # t in [1.05,300], zf in [0.02,0.98], A in [1e-3,1e6], Bf in [1,100]
                                # (many hours on a laptop; the zf > 0.9 tiles are the slow ones)
Requires: pip install mpmath
