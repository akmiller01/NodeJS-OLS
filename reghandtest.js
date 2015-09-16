var sylvester = require('sylvester'),
    ols = require('./ols.js');
    
function randomreg(n,k){
    var Y = sylvester.Matrix.Random(n,1);
    Y = Y.x(10);
    Y = Y.round();
    var X = sylvester.Matrix.Random(n,k);
    X = X.x(10);
    X = X.round();
    console.log('Y: ' + Y.elements);
    console.log('X: ' + JSON.stringify(X.transpose().elements));
    var reg = ols.reg(Y,X,true);
    return reg;
}

function staticreg() {
    var Y = $M([10,4,3,5,9,3,3,2,6,7]);
    var X = $M([[6,4,9,1,4,6,4,0,9,10],[1,4,2,5,7,3,2,9,6,1],[6,0,2,3,6,2,8,8,5,7],[2,6,8,10,1,2,8,2,6,1]]);
    X = X.transpose();

    var reg = ols.reg(Y, X);
    return reg;
}

console.log('Random: ' + JSON.stringify(randomreg(10, 4), null, 2));
// console.log('Known: ' + JSON.stringify(staticreg(), null, 2));