
biocLite("pdInfoBuilder")
naeData = getAE(
  "E-GEOD-73374",
  path = current_path,
  sourcedir=current_path,
  local = TRUE,
  type = 'raw')


install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/pd.hugene20st.hs.entrezg_22.0.0.tar.gz", repos = NULL)

oligoData = oligo::read.celfiles(filenames = filePaths, pkgname = "pd.hugene20st.hs.entrezg")
