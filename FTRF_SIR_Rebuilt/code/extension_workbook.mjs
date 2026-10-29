import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';
const root=path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const blocks=JSON.parse(await fs.readFile(path.join(root,'results/extension/workbook_data.json'),'utf8'));
const wb=Workbook.create();
const preview=path.join(root,'audit/workbook_previews'); await fs.mkdir(preview,{recursive:true});
function col(n){let s='';for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
for(const b of blocks){
  const sh=wb.worksheets.add(b.name); const end=col(b.columns.length-1); const last=b.rows.length+5;
  sh.showGridLines=false; sh.freezePanes.freezeRows(5);
  sh.getRange(`A2:${end}${last}`).format.font={name:'Arial',size:10};
  sh.getRange('A2').values=[[b.title]];
  sh.getRange('A2').format.font={name:'Arial',size:14,bold:true};
  sh.getRange('A3').values=[['Source: '+b.source]];
  sh.getRange('A4').values=[[b.note]];
  sh.getRange(`A5:${end}${last}`).values=[b.columns.map(x=>x.replaceAll('_',' ')),...b.rows];
  const body=sh.getRange(`A5:${end}${last}`); body.format.font={name:'Arial',size:10};
  body.format.columnWidth=17;body.format.rowHeight=19;
  sh.getRange(`A5:${end}5`).format={fill:'#344f6c',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:44};
  for(let j=0;j<b.columns.length;j++){
    const key=b.columns[j],c=col(j),range=sh.getRange(`${c}6:${c}${last}`);
    if(['method','reference'].includes(key))sh.getRange(`${c}5:${c}${last}`).format.columnWidth=26;
    if(['dataset','region','interval','contrast','scenario','prior','omitted_dimension','omitted_value'].includes(key))range.setNumberFormat('@');
    else if(['n','budget','horizon','season','clusters','rep'].includes(key))range.setNumberFormat('0');
    else if(['coverage95','coverage_error95','fraction_better'].includes(key))range.setNumberFormat('0.0%');
    else range.setNumberFormat('0.000');
    if(key==='rmse'&&b.columns.includes('mse')){
      const source=col(b.columns.indexOf('mse'));
      range.formulas=b.rows.map((_,i)=>[`=SQRT(${source}${i+6})`]);
    }
    if(key==='coverage_error95'&&b.columns.includes('coverage95')){
      const source=col(b.columns.indexOf('coverage95'));
      range.formulas=b.rows.map((_,i)=>[`=ABS(${source}${i+6}-0.95)`]);
    }
  }
}
wb.recalculate();
for(const b of blocks){
  const image=await wb.render({sheetName:b.name,range:`A1:${col(Math.min(b.columns.length,8)-1)}14`,scale:1.5,format:'png'});
  await fs.writeFile(path.join(preview,b.name.replaceAll(' ','_')+'.png'),new Uint8Array(await image.arrayBuffer()));
}
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!',options:{useRegex:true,maxResults:20},maxChars:1000})).ndjson);
const xlsx=await SpreadsheetFile.exportXlsx(wb); await xlsx.save(path.join(root,'results/extension_results.xlsx'));
console.log('Workbook exported, 6 sheets');
