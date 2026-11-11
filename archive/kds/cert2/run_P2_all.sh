#!/bin/sh
# Item 121: exact 3-D Bernstein certificate for the (P2) lemma, charts covering [0,1]^3 in (T = 1/t, zeta, sigma = S/x0)
for ch in R1m R1n R1p R1q R2p R2q R2m_x R2m_y R2m_p R2m_q R2n_x R2n_y R2n_p R2n_q; do python3 p2cert.py $ch; done
for ch in A B_x B_y B_p B_q; do python3 p2cert.py $ch 1/10; done
