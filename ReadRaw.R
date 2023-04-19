library(rawrr)
#basePeak gives you base peak intensity of a spectrum
rawFile = "Y:\\2023\\4\\Greedy Data\\LF_HeLa_HCD23_BE2_rep1_CH2.raw"
x = readSpectrum(
rawfile = rawFile,
scan = c(60652),
tmpdir = tempdir(),
validate = FALSE,
mode = ""
)
print(basePeak(x[[1]])[2])
