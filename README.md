# Roth Tensor Codes
Here is a collection of data and MAGMA code associated with the study of the generalization of the tensor codes introduced by Ron Roth. 

## Table of contents
:open_file_folder: <ins>__Counting elements__</ins>\
**ComputingCardinalities.txt**: *MAGMA code to compute the number of words with the conditions on the rank, the tensor-rank ect..*.


:open_file_folder: <ins>__Decoding algorithms__</ins>\
**FibreWiseDecoders.m**: *MAGMA code to decode a tensor codeword with the fibrewise decoder*.\
The algorithm ```ColumnWiseDecoder(n,mu1,mu2)```, i.e. the fibre-wise decoder only in one direction.\
The algorithm ```FibreWiseDecoder(n,mu1,mu2)```, i.e. the fibre-wise decoder with correction in the second direction.\
The algorithm used for Gabidulin decoding is in [this paper by Antonia Wachter-Zeh, Valentin Afanassiev & Vladimir Sidorenko](https://link.springer.com/content/pdf/10.1007/s10623-012-9659-5.pdf?pdf=inline%20link).\
     - 
