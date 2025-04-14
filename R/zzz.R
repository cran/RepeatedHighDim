.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "
  __        __   _                            _
  \\ \\      / /__| | ___ ___  _ __ ___   ___  | |_ ___
   \\ \\ /\\ / / _ \\ |/ __/ _ \\| '_ ` _ \\ / _ \\ | __/ _ \\
    \\ V  V /  __/ | (_| (_) | | | | | |  __/ | || (_) |
     \\_/\\_/ \\___|_|\\___\\___/|_| |_| |_|\\___|  \\__\\___/

    ",
  "     ____                       _           _
        |  _ \\ ___ _ __   ___  __ _| |_ ___  __| |
        | |_) / _ \\ '_ \\ / _ \\/ _` | __/ _ \\/ _` |
        |  _ <  __/ |_) |  __/ (_| | ||  __/ (_| |
        |_| \\_\\___| .__/ \\___|\\__,_|\\__\\___|\\__,_|
                  |_| _   _ _       _
                     | | | (_) __ _| |__
                     | |_| | |/ _` | '_ \\
                     |  _  | | (_| | | | |
                     |_| |_|_|\\__, |_| |_|  _
                        ____  _|___/       | |
                       |  _ \\(_)_ __ ___   | |
                       | | | | | '_ ` _ \\  |_|
                       | |_| | | | | | | |  _
                       |____/|_|_| |_| |_| (_)

    \n",

    "\n",
    "Toolkit for analyzing high-dimensional repeated measurements,\n",
    "providing functions for:\n",
    "- Binary random data generation\n",
    "- Differential expression analysis\n",
    "- Gene-set tests\n",
    "- Network meta-analysis for gene expression data\n",
    "- Outlier detection\n",
    "- ScRNA-seq clustering with outlier detection\n"
  )

  packageStartupMessage("This is version ", utils::packageVersion(pkgname),
                        " of ", pkgname,".")
}
