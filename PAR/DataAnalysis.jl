using NCDatasets
nc = Dataset("./PAR/Data/AQUA_MODIS.20230101.L3b.DAY.PAR.x.nc")
par = nc["BinList"]
