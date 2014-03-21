var sylvester = require('sylvester'),
    ols = require('./ols.js');
    
function randomreg(n,k){
    var Y = sylvester.Matrix.Random(n,1),
    X = sylvester.Matrix.Random(n,k),
    time = ols.reg(Y,X,true).overall.Time;
    return time;
}

console.log("Observation test...");
var sec = 0,
    i = 0,
    step = 50,
    maxobs = 0-step,
    obslimit = 1;
while (sec<=obslimit){
    maxobs+= step;
    i+=1;
    var n = (i)*step,
    k = 1,
    sec = randomreg(n,k);
    console.log("N: "+n+", K:"+k+", Sec: "+sec);
}

console.log("Param test...");
var sec = 0,
    i = 0,
    step = 5,
    maxparams = 0-step,
    paramslimit = 1;
while (sec<=paramslimit){
    maxparams+= step;
    i+=1;
    var n = (i*step)+2,
    k = i*step,
    sec = randomreg(n,k);
    console.log("N: "+n+", K:"+k+", Sec: "+sec);
}

console.log("Maximum observations in "+obslimit+" second(s) with one parameter: " + maxobs);
console.log("Maximum parameters in "+obslimit+" second(s) with minimum observations: " + maxparams);
console.log("Stress test complete.");