# Overview
Metajet is a header-only C++ library providing a fast O(N^2) template-based
implementation of generalized kt jet clustering algorithms [[1][ref1],[2][ref2]].
In contrast to FastJet [[3][ref3]], metajet does not provide O(NlogN) algorithms.
Our goal is to provide a more flexible, more generally applicable, and more
optimized implementation of O(N^2) algorithms.

[ref1]: http://arxiv.org/abs/0802.1189
[ref2]: http://arxiv.org/abs/hep-ph/0512210
[ref3]: http://arxiv.org/abs/1111.6097

# Installation
Copy the [`metajet.hh`](/src/metajet.hh) file into your project. That's it!

# Remarks
1. If you are not compiling with GCC, you may have to get rid of
   `__builtin_expect()`s. On a unix system, you can do this with a command like
   `sed -i 's/__builtin_expect(\(.*\),[01])/\1/g' metajet.hh`.
   If your compiler provides a similar facility, you may want to substitute
   the respective statement instead.
