# AlphaCore

These are the supplementary files for the AlphaCore KDD paper.

In order to run, you first need to download the data (roughly 3GB), hosted at: https://mega.nz/file/6WIiEQ5b#s4HmjOkO9WbB_-egCpmjqWoRup9iVUG2RaR6IBtS3UA
Next, decompress it, and store it in the folder data/tokens

Then open up the file evaluation.R from the main directory.

This is the main file to run the evaluation.
Input your number of CPU cores and the path to the extracted database file.

To run the entire evaluation, ideally do so on a server with Rscript, as some algorithms, i.e. some centralities depending on all pairs shortest path computations, and our own implementation of weighted k-core have very long running times. The full executing takes between 1-3 days, depending on hardware.

AlphaCore with an exponentially decaying step size however is quite fast, even though it is only an R implementation.
