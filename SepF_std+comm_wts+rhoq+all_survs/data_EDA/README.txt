catch.dat - data sent by DFO
RV_IO.dat - DFO research Vessel (RV) mean number per tow for offshore+inshore strata. Data sent by DFO
RV_O.dat - DFO research Vessel (RV) mean number per tow for offshore strata. Data sent by DFO
Landings.dat - extracted from DFO Can. Sci. Advis. Sec. Res. Doc. 2017/063. Data prior to 1980 extract from an earlier report.
mat.dat - extracted from DFO Can. Sci. Advis. Sec. Res. Doc. 2017/063.

midy_wt.dat and stock_wt.dat - output from a cohort model but I updated length~weight relationship based on RV data in DFO Can. Sci. Advis. Sec. Res. Doc. 2017/063
** there seems to be some between-year variation in length-weight relationship but I did not pursue this.
** Some evidence that a has decreased over years, and b increased, where W=aL^b. A model with AR(1) deviations in a and b may be useful.

comm_wt.txt - copied from DFO Can. Sci. Advis. Sec. Res. Doc. 2017/063.
DFO_stock_wt.txt - data sent by DFO.
These last 2 files seem unreliable to me and I suggest you dont use them directly - but catch weights need more thinking.
Also, somehow you need to extrapolate stock and catch weights from 1982 back to 1959.