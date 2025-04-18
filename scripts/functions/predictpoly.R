
# function to get predicted values from a polynomial or spline fit, recentered where desired

predictpoly <- function(model,xval,fit="poly",order=2,knots=NA,df=NA,recenter=NA) {
  if (fit=="poly") {
    pred = as.numeric(t(as.matrix(coef(model)[1:order]))%*%t(matrix(nrow=length(xval),ncol=order,data=poly(xval,order,raw=T))))
  }
  if (fit=="spline") {
    if (is.na(knots)==F) {
      ll=length(knots)+1
      pred = as.numeric(t(as.matrix(coef(model)[1:ll]))%*%t(matrix(nrow=length(xval),ncol=ll,data=ns(xx,knots=knots))))
    } else {
      pred=as.numeric(t(as.matrix(coef(model)[1:df]))%*%t(matrix(nrow=length(xval),ncol=df,data=ns(xx,df=df))))
    }
  }
  if (is.na(recenter)==F) {
    pred = pred - pred[xval==recenter]
  }
}
