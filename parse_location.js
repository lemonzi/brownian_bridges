#!/usr/bin/env node

var data = require('./locations.json');
data.locations.forEach(function(l) {
   console.log([l.timestampMs, l.latitudeE7, l.longitudeE7, l.accuracy].join('\t'));
});

