var sylvester = require('sylvester');
var Pi=Math.PI,
    PiD2=Pi/2;
function StatCom(q,i,j,b) {
    var zz=1; var z=zz; var k=i; while(k<=j) { zz=zz*q*k/(k-b); z=z+zz; k=k+2;}
    return z;
    }
function StudT(t,n) {
    t=Math.abs(t); var w=t/Math.sqrt(n); var th=Math.atan(w);
    if(n==1) { return 1-th/PiD2;}
    var sth=Math.sin(th); var cth=Math.cos(th);
    if((n%2)==1)
            { return 1-(th+sth*cth*StatCom(cth*cth,2,n-3,-1))/PiD2;}
            else
            { return 1-sth*StatCom(cth*cth,1,n-3,-1);}
    }
function reg(Y, X1, Robust) {
    var startTime = new Date().getTime();
    if (Robust === undefined) Robust = false;
    var Ones = sylvester.Matrix.Ones(X1.rows(),1),
        X = Ones.augment(X1),
        n = X.rows(),
        k = X.cols();
    if ((X.transpose().x(X)).inverse() === null) return "Collinearity error";
    if (n-k<=0) return "Too few degrees of freedom for estimating unknowns";
    var XtransXinv = (X.transpose().x(X)).inverse(),
        B = XtransXinv.x((X.transpose().x(Y))),
        Yhat = X.x(B),
        E = Y.subtract(Yhat),
        S2 = (E.transpose().x(E)).x(1/(n-B.rows())).e(1,1),
        RMSE = Math.sqrt(S2),
        VarCov = XtransXinv.map(function(d){return d*S2;}),
        SE = VarCov.diagonal().map(function(d){return Math.sqrt(d);}),
        Tstat = B.col(1).toDiagonalMatrix().x(SE.toDiagonalMatrix().inverse()).diagonal(),
        Ybar = Y.col(1).toDiagonalMatrix().trace()/n,
        R2 = 1-((E.transpose().x(E)).e(1,1)/((Y.map(function(d){return d-Ybar;})).transpose().x((Y.map(function(d){return d-Ybar;})))).e(1,1)),
        E2 = E.x(E.transpose()),
        Vhat = E2.diagonal().toDiagonalMatrix(),
        RVarCov = XtransXinv.x(X.transpose()).x(Vhat.transpose()).x(X).x(XtransXinv).map(function(d){return d*(n/(n-k));}),
        RSE = RVarCov.diagonal().map(function(d){return Math.sqrt(d);}),
        RTstat = B.col(1).toDiagonalMatrix().x(RSE.toDiagonalMatrix().inverse()).diagonal(),
        result = {};
        result.overall = {'obs':n,'params':k,'R2':R2,'RMSE':RMSE,'Robust':Robust,'Time':((new Date().getTime())-startTime)/1000};
    var i = 0;
    for(i=0;i<k;i++){
        if (Robust === false){
            result['B'+i] = {'value':B.e(i+1,1),'SE':SE.e(i+1,1),'Tstat':Tstat.e(i+1,1),'Pval':StudT(Tstat.e(i+1,1),n-k)};
        } else {
            result['B'+i] = {'value':B.e(i+1,1),'SE':RSE.e(i+1,1),'Tstat':RTstat.e(i+1,1),'Pval':StudT(RTstat.e(i+1,1),n-k)};
        }
    }
    return result;
}

exports.reg = reg;