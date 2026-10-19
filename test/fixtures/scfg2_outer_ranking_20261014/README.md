# Outer ranking integration witness

`bpRNA_RFAM_27342.txt` is the 132-nt given-scaffold input already used in the
2026-10-14 pilot-exponent and constant-factor studies. It is frozen here so
API/CLI tests do not depend on ignored benchmark artifacts.

At baseline `9c76f292`, Turner04 / dangles=2 / RequireWeakClosure / beam20:
ordinary ranking gives log Z=56.550132046364084; outer with q=2/3, p=8 and
frontier4p gives log Z=60.715913323849691. These are ranking-regression values,
not an independent semantic oracle. Patched legacy inside remains primary.
