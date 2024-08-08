#### download data from Zenodo ####
# download latest version of the database from Zenodo, https://doi.org/10.5281/zenodo.7494930

cat("Downloading raw data\n")
zen4R::download_zenodo("10.5281/zenodo.13270547", path = "data/smlbase/")

cat("Data successfully downloaded! \n")
