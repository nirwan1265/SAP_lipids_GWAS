#https://github.com/Waddlessss/MAFFIN/

devtools::install_github("Waddlessss/MAFFIN")


library(MAFFIN)

View(TestingData)
inputTesting = TestingData

# Sample normalization using MAFFIN main function
MAFFINTable = MAFFINNorm(inputTesting)
