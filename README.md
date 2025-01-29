# Roth Tensor Codes
Here is a collection of data and MAGMA code associated with the study of the generalization of the tensor codes introduced by [Ron M. Roth](https://ieeexplore.ieee.org/document/556603). One can find on [this official MAGMA webpage](https://magma.maths.usyd.edu.au/calc/) a free magma calculator. All the functions in the code, except if trivial or repeated, have a linked procedure to test them.

## Table of contents
:open_file_folder: <ins>__Counting elements__</ins>\
**ComputingCardinalities.txt**: *MAGMA code to compute the number of 3-tensor of tensor-rank 1 or 2, as well as the number of matrices of a given rank, and the number of errors corrected by the fibre-wised decoders below*.

:open_file_folder: <ins>__Decoding algorithms__</ins>\
**FibreWiseDecoders.m**: *MAGMA code to decode a tensor codeword with the fibrewise decoders*.\
The algorithm ```ColumnWiseDecoder(R,mu1,mu2)```, i.e. the fibre-wise decoder only in one direction.\
The algorithm ```FibreWiseDecoder(R,mu1,mu2)```, i.e. the fibre-wise decoder with correction in the second direction.\
The algorithm used for Gabidulin decoding is in [this paper by Antonia Wachter-Zeh, Valentin Afanassiev & Vladimir Sidorenko](https://link.springer.com/content/pdf/10.1007/s10623-012-9659-5.pdf?pdf=inline%20link).\
\
**RadicalDecoders.m**: *MAGMA code to decode a tensor codeword with the radical decoders*.\
The algorithm ```FactoringOnTheLeft(N,V,S)``` factors any q-polynomial N(X,Y) = V(f(X,Y)) with V(Z) a single variable q-polynomial.\
The algorithm ```RadicalDecoderFix(R,mu,t)```, decodes a codeword with an error thas has 3-fibre weight upper bounded by t using the radical criterion for factorization.\
The algorithm ```RadicalDecoder(R,mu)```, decodes a codeword with an error that has wSigma upper bounded, using the radical criterion for factorization.

**DecodingInTensorRankMetrix.m**: *MAGMA code to decode a tensor codeword with an error with tensor-rank upper bounded*.\
The algorithm ```RandomErrorTRANK(T)``` creates a random (non-uniform) matrix over Fqn whith tensor rank (over Fq) upper-bounded by T.\
The algorithm ```DecodeTRANKbyProxy(R,mu)``` decodes a codeword with an error whose tensor-rank is upper-bounded.\
\
\
:hammer_and_wrench: <ins>*To be added:*</ins> *Generalizations to order m-tensors, syndrome tensor-rank decoder upgraded.*
