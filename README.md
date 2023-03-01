# DATESModifiedCurveFit
Find 2 admixture dates from DATES output

Run modified fits on various DATES output folders
Description: Estimate admixture timing for 1 or 2 admixtures

Usage:

runfits(
  pvalue = 0.001,
  directorylist = NULL,
  foldernameformat = "NoFormat",
  outputfilename = NULL
)

Arguments

pvalue: P-value below which the program attempts a double exponential decay fit. Default 0.001

directorylist: List of directories with DATES output. Use list.dirs function on R to generate (remove the first parent folder from list). Ensure all folders have DATES output. Required field.

foldernameformat: NoFormat = any folder name of choice, any other input = T1_T2_adm1_adm2 format for folder name where T1,T2 are admixture times (from simulation) and adm1,2 are admixture percentages in each wave (from simulation). Default = "NoFormat"

outputfilename: Optional, Specify file name. Output is written to file in TSV format.

Output: A dataframe with DATES exponential fitted output for 1 admixture (and 2 in case of p-value below threshold)

Examples:

out1 <- runfits(pvalue=0.001,directorylist=dirs,foldernameformat="NoFormat");
out1 <- runfits(directorylist=dirs,foldernameformat="Any",outputfilename="/home/user/output.tsv");
