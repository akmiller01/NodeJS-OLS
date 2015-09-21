var sylvester = require('sylvester');
var Pi = Math.PI,
    PiD2 = Pi / 2;

function StatCom(q, i, j, b) {
    var zz = 1;
    var z = zz;
    var k = i;
    while (k <= j) {
        zz = zz * q * k / (k - b);
        z = z + zz;
        k = k + 2;
    }
    return z;
}

function StudT(t, n) {
    t = Math.abs(t);
    var w = t / Math.sqrt(n);
    var th = Math.atan(w);
    if (n == 1) {
        return 1 - th / PiD2;
    }
    var sth = Math.sin(th);
    var cth = Math.cos(th);
    if ((n % 2) == 1) {
        return 1 - (th + sth * cth * StatCom(cth * cth, 2, n - 3, -1)) / PiD2;
    } else {
        return 1 - sth * StatCom(cth * cth, 1, n - 3, -1);
    }
}

function FishF(f, n1, n2) {
    var x = n2 / (n1 * f + n2)
    if ((n1 % 2) == 0) {
        return StatCom(1 - x, n2, n1 + n2 - 4, n2 - 2) * Math.pow(x, n2 / 2)
    }
    if ((n2 % 2) == 0) {
        return 1 - StatCom(x, n1, n1 + n2 - 4, n1 - 2) * Math.pow(1 - x, n1 / 2)
    }
    var th = Math.atan(Math.sqrt(n1 * f / n2));
    var a = th / PiD2;
    var sth = Math.sin(th);
    var cth = Math.cos(th)
    if (n2 > 1) {
        a = a + sth * cth * StatCom(cth * cth, 2, n2 - 3, -1) / PiD2
    }
    if (n1 == 1) {
        return 1 - a
    }
    var c = 4 * StatCom(sth * sth, n2 + 1, n1 + n2 - 4, n2 - 2) * sth * Math.pow(cth, n2) / Pi
    if (n2 == 1) {
        return 1 - a + c / 2
    }
    var k = 2;
    while (k <= (n2 - 1) / 2) {
        c = c * k / (k - .5);
        k = k + 1
    }
    return 1 - a + c
}

function StatCom(q, i, j, b) {
    var zz = 1;
    var z = zz;
    var k = i;
    while (k <= j) {
        zz = zz * q * k / (k - b);
        z = z + zz;
        k = k + 2
    }
    return z
}

function reg(Y, X1, parameters) {
    var startTime = new Date().getTime();
    var Ones = sylvester.Matrix.Ones(X1.rows(), 1);
    var X = Ones.augment(X1);
    var n = X.rows();
    var k = X.cols();
    var XtransXinv = (X.transpose().x(X)).inverse();
    if (XtransXinv === null) return "Collinearity error";
    if (k >= n) return "Too few degrees of freedom for estimating unknowns (" + k + " variables but only " + n + " observations)";
    var B = XtransXinv.x((X.transpose().x(Y)));
    var Yhat = X.x(B);
    var E = Y.subtract(Yhat);
    var S2 = (E.transpose().x(E)).x(1 / (n - B.rows())).e(1, 1);
    var RMSE = Math.sqrt(S2);
    var VarCov = XtransXinv.map(function(d) {
        return d * S2;
    });
    var SE = VarCov.diagonal().map(function(d) {
        return Math.sqrt(d);
    });
    var diagonalSEInverse = SE.toDiagonalMatrix().inverse();
    if (diagonalSEInverse == null) {
        return "Could not calculate Tstat"
    }
    var Tstat = B.col(1).toDiagonalMatrix().x(SE.toDiagonalMatrix().inverse()).diagonal();
    var Ybar = Y.col(1).toDiagonalMatrix().trace() / n;
    var SST = ((Y.map(function(d) {
        return d - Ybar;
    })).transpose().x((Y.map(function(d) {
        return d - Ybar;
    })))).e(1, 1);
    var R2 = 1 - ((E.transpose().x(E)).e(1, 1) / SST);
    var Fstat = ((R2 * SST) / (k - 1)) / S2;

    var result = {};
    result.overall = {
        'obs': n,
        'params': k,
        'R2': R2,
        'AdjR2': ((1 - R2)*(n - 1)) / (n - k),
        'RMSE': RMSE,
        'Fstat': Fstat,
        'Fvalue': FishF(Fstat, k - 1, n - k),
        'Time': ((new Date().getTime()) - startTime) / 1000
    };
    var i = 0;
    for (i = 0; i < k; i++) {
        var name = (parameters && i !== 0) ? parameters[i - 1] : 'B' + i;
        result[name] = {
            'value': B.e(i + 1, 1),
            'SE': SE.e(i + 1, 1),
            'Tstat': Tstat.e(i + 1, 1),
            'Pval': StudT(Tstat.e(i + 1, 1), n - k)
        };
    }
    return result;
}

exports.reg = reg;