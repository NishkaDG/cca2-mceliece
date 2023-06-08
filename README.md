# cca2-mceliece
A complete implementation of the McEliece system in its original form and CCA2-secure versions

These files implement the Classic McEliece PKC as well as conversions of it to IND-CCA2 security using the Fujisaki-Okamoto transform and Kobara-Imai alpha transform. 
Along the way it implements algorithms to convert binary strings to vectors of length n and weight t (also called error vectors for their use in (n,t)-error correcting codes).

## Requirements
Python3
[SageMath](https://www.sagemath.org/)
[Pycryptodome](https://pycryptodome.readthedocs.io/en/latest/)

## Detailed Description
For a general overview of McEliece-based cryptosystems, click [here](http://classic.mceliece.org/)
### classic.py
This file implements [Classic McEliece](https://ipnpr.jpl.nasa.gov/progress_report2/42-44/44N.PDF). 
Decoding of the underlying Goppa code is done with code in bernstein.py, sourced from [https://cr.yp.to/papers/goppadecoding-20220816.pdf](https://cr.yp.to/papers/goppadecoding-20220816.pdf) (auto-converted from Sage to Python).

### cca_conversions.py
This file implements the [Fujisaki-Okamoto transform](https://link.springer.com/content/pdf/10.1007/s00145-011-9114-1.pdf), Cayrel et al's [variant](https://hal-ujm.archives-ouvertes.fr/file/index/docid/712875/filename/2012_PKC_cayrel.pdf) on the Fujisaki-Okamoto transform, and [Kobara-Imai alpha transform](https://link.springer.com/content/pdf/10.1007/3-540-44586-2_2.pdf).
Cayrel et al use Srivastava codes in their proposal but I use Goppa codes in my implementation for stronger security guarantees.
The other CCA conversions require PRNGs and functions to convert binary strings to error vectors; these are implemented in other files.

### auxiliary.py 
This file implements the hash functions/PRNGs using cshake in pycryptodome, as well as other helper functions.

### sendrier.py
This file implements [Sendrier's protocol](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1523371&tag=1) for converting binary strings into error vectors for given parameters of n, t.

### ideal_stc.py
This file implements [Barenghi and Pelosi's protocol](https://re.public.polimi.it/bitstream/11311/1137353/3/3387902.3392630.pdf) for converting binary strings into error vectors for given parameters of n, t.

## Bibliography 
(In case the above links break)
Classic McEliece: Robert J. McEliece. “A Public-Key Cryptosystem Based On Algebraic Coding Theory”. In: Deep Space Network Progress Report 44 (Jan. 1978), pp. 114–116. url: https://ipnpr.jpl.nasa.gov/progress_report2/42-44/44N.PDF
The Fujisaki-Okamoto generic transform: Eiichiro Fujisaki and Tatsuaki Okamoto. “Secure integration of asymmetric and symmetric encryption schemes”. In: Advances in Cryptology—CRYPTO’99: 19th Annual International Cryptology Conference Santa Barbara, California, USA, August 15–19, 1999 Proceedings. Springer. 1999, pp. 537–554.
Cayrel et al's variant of the Fujisaki-Okamoto transform: Pierre-Louis Cayrel, Gerhard Hoffmann, and Edoardo Persichetti. “Efficient Implementation of a CCA2-Secure Variant of McEliece Using Generalized Srivastava Codes”. In: International Conference on Theory and Practice of Public Key Cryptography. 2012.
Kobara-Imai: Kazukuni Kobara and Hideki Imai. “Semantically Secure McEliece Public-Key Cryptosystems - Conversions for McEliece PKC”. In: International Conference on Theory and Practice of Public Key Cryptography. 2001
Bernstein's Goppa decoding paper: Daniel J. Bernstein. Understanding binary-Goppa decoding. Cryptology ePrint Archive,Paper 2022/473. url: https://eprint.iacr.org/2022/473. 2022.
Sendrier's conversion: Nicolas Sendrier. “Encoding information into constant weight words”. In: Proceedings. International Symposium on Information Theory, 2005. ISIT 2005. 2005, pp. 435–438. doi: 10.1109/ISIT.2005.1523371
Barenghi and Pelosi's conversion: Alessandro Barenghi and Gerardo Pelosi. “Constant weight strings in constant time: a building block for code-based post-quantum cryptosystems”. In: Proceedings of the 17th ACM International Conference on Computing Frontiers (2020)

See also: D. Engelbert, R. Overbeck, and A. Schmidt. A Summary of McEliece-Type Cryptosystems and their Security. Cryptology ePrint Archive, Paper 2006/162. url: https://eprint.iacr.org/2006/162.
