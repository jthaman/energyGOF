# Method for adding a new distribution

1. Write the new S3 GOFDist object.
2. Write new Roxygen block.
3. Add line to table in -package.R
4. Update DESCRIPTION
5. Update formals of energyGOF.test
6. Update char_to_dist function.
7. Update readme to list new dist
8. Write unit tests.

# Goals to submit package to CRAN

1. Pass cran checks  (error with fitdistrplus)
2. (DONE) Check documentation for spelling and grammar
3. Better unit testing
4. Power analysis in a vignette
