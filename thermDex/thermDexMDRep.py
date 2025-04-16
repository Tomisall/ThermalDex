import sys
import subprocess
from pathlib import Path
#from markdown_pdf import MarkdownPdf, Section
from thermDex.custom_Markdown_PDF import MarkdownPdf, Section
import fitz
#from thermDexMolecule import thermalDexMolecule


def mdReportCreation(molecule, dataURL):
    memo = MarkdownPdf(toc_level=0)

    #csspath = Path("./_core/style.css")
    #css = csspath.read_bytes().decode()
    css = '''
table {
    border-spacing: 2px;
    width: 100%;
    table-layout: fixed;
    border-collapse: collapse;
}

td {
    padding: 8px;
    border: 1px solid;
}

th {
    padding: 8px;
    border: 1px solid;
}

tr {
    background-color: #fff;
    border-top: 1px solid;
}
'''
    if molecule.mp == '' and molecule.mpEnd == '':
        mpString = ''

    elif molecule.mp != '' and molecule.mpEnd == '':
        mpString = f'mp: {molecule.mp} °C'

    elif molecule.mp != '' and molecule.mpEnd != '':
        mpString = f'mp: {molecule.mp} to {molecule.mpEnd} °C'

    report_body = f'''
<html>
<head>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=UTF-8">
</head>
<body>
<div style="text-align: center;">
<h1>Thermal Hazard Assessment Memo</h1>
</div>
<h2>{molecule.name}</h2>
<div style="text-align: center;">
<img src="{dataURL}" alt="HTML image test" width="50%"/>
</div>              
<h3>Molecule Properties</h3>
<p>SMILES: {molecule.SMILES}</p>
<p>Formula: {molecule.eleComp}</p>
<p>{mpString}</p>
<h3>Results</h3>

<table align="Center" style="border-spacing: 10px;">
<thead>
<tr class="head"><td></td><td></td><td></td></tr>
</thead>
<tbody>
<tr class="even">
<td colspan="3">High Energy Groups: ({molecule.HEG}) {molecule.HEG_list} &nbsp;</td>
</tr>
<tr class="even">
<td colspan="3">Explosive Groups: ({molecule.EFG}) {molecule.EFG_list}</td>
</tr>
<tr class="even">
<td>Rule of Six = {molecule.RoS_val}</td>
<td>Oxygen Balance = {molecule.OB_val}</td>
<td> </td>
</tr>
<tr class="even">
<td>Q<sub>DSC</sub> = {molecule.Q_dsc} {molecule.Qunits}</td>
<td>T<sub>onset</sub> = {molecule.onsetT}</td>
<td>T<sub>init</sub> = {molecule.initT}</td>
</tr>
<tr class="even">
<td>Impact Sensitivity = {molecule.IS_val}</td>
<td>Explosive Propagation = {molecule.EP_val}</td>
<td>T<sub>D24</sub> = {molecule.Td24} °C</td>
</tr>
</tbody>
</table>
<p>O.R.E.O.S. assessment of risk by scale:</p>
<table>
<thead>
<tr class="header">
<th> <5 g</th>
<th> 5 to 100 g</th>
<th> 100 to 500 g</th>
<th> >500 g</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>{molecule.oreoSmallScale_des}</td>
<td>{molecule.oreoTensScale_des}</td>
<td>{molecule.oreoHundsScale_des}</td>
<td>{molecule.oreoLargeScale_des}</td>
</tr>
</tbody>
</table>
<h3>Interpretation</h3>
<p>These results have been calculated using X<sup>1</sup> and they show Y<sup>2</sup>.</p>
<p><small>[1]: <i>Org. Proc. Res. Dev.,</i> 2011, 2341-2356</small><br>
<small>[2]: <i>Org. Proc. Res. Dev.,</i> 2011, 2117-2119</small></p>
</body>
'''

    memo.add_section(Section(report_body, toc=False), css)

    memo.meta.update({
      "creationDate": fitz.get_pdf_now(),
      "modDate": fitz.get_pdf_now(),
      "creator": "ThermalDex: Thermal Hazard Assessment Tool",
      "producer": None,
      "title": "Thermal Hazard Assessment Memo",
      "author": "ThermalDex: Thermal Hazard Assessment Tool",
      "subject": None,
      "keywords": None,
    })

    memo.save("./AssessmentMemos/MD_thermal_assessment_memo.pdf")

    doc = fitz.open("./AssessmentMemos/MD_thermal_assessment_memo.pdf")
    numpages = doc.page_count  # number of pages
    footer_text = "Page %i of %i"
    header_text = "My Matplotlib File"
    blue = fitz.pdfcolor["blue"]

    for page in doc:
        prect = page.rect
        header_rect = fitz.Rect(0, 36, prect.width, 56)  # height 20 points
        page.insert_textbox(header_rect, header_text,
                            fontname="hebo", color=blue,
                            align=fitz.TEXT_ALIGN_CENTER)

        ftext = footer_text % (page.number + 1, numpages)
        y1 = prect.height - 36  # bottom of footer rect
        y0 = y1 - 20  # top of footer rect
        footer_rect = fitz.Rect(0, y0, prect.width, y1)  # rect has full page width
        page.insert_textbox(footer_rect, 'This report may contain confidential information', align=fitz.TEXT_ALIGN_CENTER)

    doc.save("./AssessmentMemos/ALT_MD_thermal_assessment_memo_ALT.pdf")

    #cwd = "./"
    #filepath =  cwd + filename #(cwd / filename).resolve()
    #print(filepath)
    subprocess.run(['start', "./AssessmentMemos/MD_thermal_assessment_memo.pdf"], check=True, shell=True)

#molecule = thermalDexMolecule(SMILES='c1ccc(C)cc1',name='Test mol Toluene')
#molecule.genMol()
#molecule.molToIMG() #.genAllValues() #highEnergyGroups,expEnergyGroups)
#imageData = molecule.molToBytes()
#dataURL = 'data:image/png;base64,' + imageData
#print(dataURL)
#mdReportCreation(molecule, dataURL)
