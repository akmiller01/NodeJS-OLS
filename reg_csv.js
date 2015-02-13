var fs = require('fs'),
    sylvester = require('sylvester'),
    prompt = require('prompt'),
    csv = require('fast-csv'),
    ols = require('./ols.js');
    
var stream = fs.createReadStream(process.argv[2]),
data = [];
 
csv.fromStream(stream, {headers : true})
.on("data", function(d){
    data.push(d);
})
.on("end", function(){
    pickVars(data)
});

function pickVars(data){
    var schema = {
        properties: {
          Y: {
            description: "Select Y variable using index",
            type:'number',
            pattern: /^[0-9]+$/,
            message: 'Please use the index number to select your Y variable',
            required: true
          },
          X: {
            description: "Select X variables using index. CTRL+C to stop",
            type:'array',
            minItems:1,
            pattern: /^[0-9]+$/,
            message: 'Please use the index numbers to select at least one X variable',
            required: true
          },
          robust: {
            description: "Robust standard errors (true or false)",
            type:'boolean',
            message: 'true or false',
            required: true
          },
        }
    };
    var headers = Object.keys(data[0]),
    headerLen = headers.length,
    dataLen = data.length;
    for(var i = 0; i < headerLen; i++){
        console.log(i+". "+headers[i]);
    };
    prompt.start();
    prompt.get(schema,function(err,result){
        var yKey = headers[result.Y],
        xKeys = [],
        xLen = result.X.length;
        for(var i = 0; i < xLen; i++){
            var thisXkey = headers[parseInt(result.X[i])];
            if(thisXkey){
                xKeys.push(thisXkey)
            };
        };
        var robust = result.robust,
        Xarr = [],
        Yarr = [];
        for(var i = 0; i < dataLen; i++){
            Yarr.push([data[i][yKey]]);
            var Xrow = [];
            for(var j = 0; j < xLen; j++){
                Xrow.push(data[i][xKeys[j]]);
            };
            Xarr.push(Xrow);
        };
        var X = $M(Xarr),
        Y = $M(Yarr);
        console.log(ols.reg(X,Y,robust));
    });
};
