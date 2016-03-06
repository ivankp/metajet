# Overview
Metajet is a header-only C++ library providing a fast $O(N^2)$ template-based
implementation of generalized $k_\mathrm{t}$ jet clustering algorithms [1,2].
In contrast to FastJet [3], metajet does not provide $O(N\log N)$ algorithms.
Our goal is to provide a more flexible, more generally applicable, and more optimized
implementation of $O(N^2)$ algorithms.

[1]: http://arxiv.org/abs/0802.1189
[2]: http://arxiv.org/abs/hep-ph/0512210
[3]: http://arxiv.org/abs/1111.6097

# Installation
Copy the `metajet.hh` file into your project. That's it!

# Remarks
1. If you are not compiling with GCC, you will have to get rid of
   `__builtin_expect()`s. On a unix system, you can do this with a command like
   ```
   sed -i 's/__builtin_expect//g' metajet.hh
   ```
