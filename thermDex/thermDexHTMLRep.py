import sys
import subprocess
from pathlib import Path
from xhtml2pdf import pisa
import time
from os import listdir, remove
from os.path import isfile, join, abspath
import fitz
from PIL import Image
from io import BytesIO
import base64
from docx2pdf import convert
from datetime import datetime

if sys.platform != "darwin":
    import win32com.client
    import pythoncom

def convert_html_to_pdf(html_string, pdf_path):
    with open(pdf_path, "wb") as pdf_file:
        pisa_status = pisa.CreatePDF(html_string, dest=pdf_file)
        
    return not pisa_status.err

def mdReportCreation(molecule, dataURL, td24Upper, td24Lower):
    wdFormatPDF = 17

    if molecule.mp == '' and molecule.mpEnd == '':
        mpString = ''

    elif molecule.mp != '' and molecule.mpEnd == '':
        mpString = f'mp: {molecule.mp} °C'

    elif molecule.mp != '' and molecule.mpEnd != '':
        mpString = f'mp: {molecule.mp} to {molecule.mpEnd} °C'


    additionalData = '''
<div class="pagebreak"> </div>
<h2>Additional Data</h2>'''
    try:
        attachedFiles = [f for f in listdir(molecule.dataFolder) if isfile(join(molecule.dataFolder, f))]
    except:
        attachedFiles = []
    print('\n\n')
    print(molecule.dataFolder)
    print(attachedFiles)
    imgList = ['.png', '.PNG', '.jpeg', '.JPEG', '.jpg', '.JPG', '.svg', '.SVG']
    pdfList = ['.pdf', '.PDF']
    docList = ['.doc', '.docx']
    for file in attachedFiles:
        if file.endswith(tuple(imgList)):
             print(file)
             additionalData += f'''
<h3>{file}</h3>
<div class="attachedFile">
<img src="{molecule.dataFolder}/{file}" alt="{file}" style="object-fit: cover;"/>
</div>
'''
        elif file.endswith(tuple(pdfList)):
            additionalData += f'<h3>{file}</h3>'
            pdfFile = fitz.open(f'{molecule.dataFolder}/{file}')
            noOfPages = pdfFile.page_count
            for pageNo in range(noOfPages):
                page = pdfFile.load_page(pageNo) 
                pix = page.get_pixmap()
                #image_data = pix.tobytes()
                #byte_array = BytesIO()
                #byte_array.write(image_data, format='PNG')
                image_bytes = pix.tobytes("PNG")
                pngOfPage = base64.b64encode(image_bytes).decode('utf-8')
                pageDataURL = 'data:image/png;base64,' + pngOfPage
                additionalData += f'''
<div class="attachedFile">
<img src="{pageDataURL}" alt="{file}" style="object-fit: cover;"/>
</div>
'''

        elif file.endswith(tuple(docList)):
            additionalData += f'<h3>{file}</h3>'
            now = datetime.now()
            neatNow = now.strftime("%d-%b-%Y_%H-%M-%S")
            outName = f'./_core/UserAddedData/temp/temp_{str(neatNow)}.pdf'
            #convert(f'{molecule.dataFolder}/{file}', outName)
            absolutePath = abspath(f'{molecule.dataFolder}/{file}')
            outAbs = abspath(outName)
            word = win32com.client.Dispatch('Word.Application', pythoncom.CoInitialize())
            doc = word.Documents.Open(absolutePath)
            doc.SaveAs(outAbs, FileFormat=wdFormatPDF)
            doc.Close()
            #word.Quit()
            pdfFile = fitz.open(outName)
            noOfPages = pdfFile.page_count
            for pageNo in range(noOfPages):
                page = pdfFile.load_page(pageNo) 
                pix = page.get_pixmap()
                #image_data = pix.tobytes()
                #byte_array = BytesIO()
                #byte_array.write(image_data, format='PNG')
                image_bytes = pix.tobytes("PNG")
                pngOfPage = base64.b64encode(image_bytes).decode('utf-8')
                pageDataURL = 'data:image/png;base64,' + pngOfPage
                additionalData += f'''
<div class="attachedFile">
<img src="{pageDataURL}" alt="{file}" style="object-fit: cover;"/>
</div>
'''
        else:
            print(f'bad {file}')

    interpRoS = f"The 'Rule of Six'<sup>1</sup> is a simple metric that checks if there are at least 6 times as many carbon atoms as there are high energy functional groups in a molecule. The result above shows the result shows the diffrence beween 6 times the number of high energy groups vs the number of carbon atoms. As such a value of &#8804; 0 implies that a compund should be safe to handel. Here the compound was found to be {molecule.RoS_des} by this metric."
    interpOB = f"The 'Oxygen Balance'<sup>1</sup> is a method for assessing the explosive properties of molecules. The vaule is based on the carbon, hydrogen and oxygen content of a molecule and the MW. Using the oxygen balance the explosive risk can be assessed as 'High (OB -120 to 80)', 'Medium (OB 80 to 160, or -240 to -120)', or 'Low' (OB <-240 or >160). Here the risk is assessed as {molecule.OB_des}."

    if molecule.Q_dsc != None and molecule.Q_dsc != '':
        QdscInfo = f"{molecule.Q_dsc} {molecule.Qunits[:-2]}<sup>-1</sup>"
    else:
        QdscInfo = "<i>Not Given</i>"

    if molecule.onsetT != None and molecule.onsetT != '':
        onsetTInfo = f"{molecule.onsetT} °C"
    else:
        onsetTInfo = "<i>Not Given</i>"

    if molecule.initT != None and molecule.initT != '':
        initTInfo = f"{molecule.initT} °C"
    else:
        initTInfo = "<i>Not Given</i>"


    blankCellCount = 0 
    if molecule.Td24 != None:
         if td24Upper >= molecule.Td24 > td24Lower:
             Td24Formated = f"<b style='color: orange;'>{'{:.1f}'.format(molecule.Td24)} °C</b>" 
             seekApproval = ''
         elif molecule.Td24 <= td24Lower:
             Td24Formated = f"<b style='color: red;'>{'{:.1f}'.format(molecule.Td24)} °C</b>" 
             seekApproval = 'Given the low T<sub>D24</sub> value for this molecule, approval must be sought internally before it is used.'
         else:
             Td24Formated = "{:.1f}".format(molecule.Td24) + " °C"
             seekApproval = ''
         TD24row = f'<td class="secretTable">T<sub>D24</sub> = {Td24Formated}</td>'
         interpTd24 = f"T<sub>D24</sub> is the temperature at which the time to the maximum rate of a runaway reaction is 24&nbsp;h.<sup>3</sup>  A reaction where the time to maximum rate is &#8805; 24&nbsp;h is highly unlikely to develop a thermal runaway. As such a reaction temperature of <b>{Td24Formated} or below</b> would be considered appropriate for this compound.{seekApproval}"
    else:
         Td24Formated = ""
         TD24row = ""
         blankCellCount += 1
         interpTd24 = "T<sub>D24</sub> (the temperature at which the time to the maximum rate of a runaway reaction is 24&nbsp;h under adiabatic conditions.<sup>3</sup>) is estmated from T<sub>init</sub>. This value would be needed (from DSC measurement) in order to calculate this useful safetey metric."

   

    if molecule.IS_val != None:
        ISvalFormated = "{:.2f}".format(molecule.IS_val)
        ISrow = f'<td class="secretTable">Impact Sensitivity = {ISvalFormated}</td>'
        if molecule.IS_val <0:
            ISinterp = "<b>not be impact sensitive</b>"
        else:
            ISinterp = "<b>be impact senstive</b>"
    else:
        ISvalFormated = ""
        ISrow = ""
        blankCellCount += 1

    if molecule.EP_val != None:
        EPvalFormated = "{:.2f}".format(molecule.EP_val)
        EProw = f'<td class="secretTableMid">Explosive Propagation = {EPvalFormated}</td>'
        if molecule.EP_val <0:
            EPinterp = "to <b>not exhibit explosive propagation</b>"
        else:
            EPinterp = "to <b>exhibit explosive propagation</b>"
    else:
        EPvalFormated = ""
        EProw = ""
        blankCellCount += 1

    fillerCells = ''
    for cell in range(blankCellCount):
        fillerCells += '<td class="secretTableMid"> &nbsp; </td>'

    if molecule.IS_val != None and molecule.EP_val != None:
        if {molecule.yoshidaMethod} == 'Pfizer':
            interpYoshida = f"The likelyhood that a compound would exhibit impact sensitive or explosive propagation can be estimated from their 'Yoshida values'. These Yoshida values have been calculated using the {molecule.yoshidaMethod} method from the Q<sub>DSC</sub> and T<sub>init</sub> data.<sup>2</sup> If the Yoshida values are <0 the compound is not expected to exhibit that property. Here the molecule is expected to {ISinterp} and {EPinterp} based on these values."
        else:
            interpYoshida = f"The likelyhood that a compound would exhibit impact sensitive or explosive propagation can be estimated from their 'Yoshida values'. These Yoshida values have been calculated using the {molecule.yoshidaMethod} method from the Q<sub>DSC</sub> and T<sub>onset</sub> data.<sup>2</sup> If the Yoshida values are <0 the compound is not expected to exhibit that property. Here the molecule is expected to {ISinterp} and {EPinterp} based on these values."
    else:
        interpYoshida =  f"The likelyhood that a compound would exhibit impact sensitive or explosive propagation can be estimated from their 'Yoshida values'. These Yoshida values can be calculated using the Yoshida method from the Q<sub>DSC</sub> and T<sub>onset</sub> data or using the more conservitive Pfizer method (recomended) from the Q<sub>DSC</sub> and T<sub>init</sub> data.<sup>2</sup> These values would need to be gathered from DSC measurement inorder to determine these important saftey metrics."

    interpOREOS = f"The O.R.E.O.S. assessment method assigns points based on the <em>O</em>xygen Balance, <em>R</em>ule of 6 value, presence of <em>E</em>xplosive functional groups, <em>O</em>nset temperature (T<sub>onset</sub>), and reaction <em>S</em>cale.<sup>1</sup> As such two diffrent point scales are used depending on if a onset temperature is provided or not. Including T<sub>onset</sub> data from DSC measurement is preferred. The table above provides the risk categories determined for this compound at 4 diffrent scales. Low Hazard results indicate reactions that can be performed using only the standard precautions. For use of compounds at scales that indicate Medium or High risk, discussion must be fist had with the process saftey team prior to any futher action. The O.R.E.O.S. value score for each scale was found to be:"

    listSMILES = [molecule.SMILES[i:i+60] for i in range(0, len(molecule.SMILES), 60)]
    #splitableSMILES = "<wbr>".join(listSMILES)
    splitableSMILES = " ".join(listSMILES)
    print(f'\n\n{splitableSMILES}\n\n')

    html_content = f'''
        <!DOCTYPE html>
        <html>
        <head>
        <META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=UTF-8">
        </head>
        <body>
        <p><link rel="stylesheet" href="./_core/style.css"></p>
        <div style="text-align: center;">
        <h1>Thermal Hazard Assessment Memo</h1>
        </div>
        <h2>{molecule.name}</h2>
        <div class="imgDiv">
        <img src="{dataURL}" alt="HTML image test" style="object-fit: cover;"/>
        </div>              
        <h3>Molecule Properties</h3>
        <p>SMILES: {splitableSMILES} <br>
        Formula: {molecule.eleComp}<br>
        MW: {"{:.2f}".format(molecule.MW)} g mol<sup>-1</sup><br>
        {mpString}</p>
        <h3>Results</h3>
        <table align="Center" class="secretTable">
        <tr class="secretTable">
        <td colspan="3" class="secretTable">High Energy Groups: ({molecule.HEG}) {', '.join(molecule.HEG_list)} &nbsp;</td>
        </tr>
        <tr class="secretTable">
        <td colspan="3" class="secretTable">Explosive Groups: ({molecule.EFG}) {', '.join(molecule.EFG_list)}</td>
        </tr>
        <tr class="secretTable" class="secretTable">
        <td class="secretTable">Rule of Six = {molecule.RoS_val}</td>
        <td class="secretTableMid">Oxygen Balance = {"{:.2f}".format(molecule.OB_val)}</td>
        <td class="secretTable"> </td>
        </tr>
        <tr class="secretTable">
        <td class="secretTable">Q<sub>DSC</sub> = {QdscInfo}</td>
        <td class="secretTableMid">T<sub>onset</sub> = {onsetTInfo}</td>
        <td class="secretTable">T<sub>init</sub> = {initTInfo}</td>
        </tr>
        <tr class="secretTable">
        {ISrow}
        {EProw}
        {TD24row}
        {fillerCells}
        </tr>
        </table>

        <p>O.R.E.O.S. assessment of risk by scale:</p>
        <table align="Center" width="65%" >
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
        <p style="margin-bottom:0%;">Emperical Test Results:</p>
        <table align="Center" class="secretTable">
        <tr class="secretTable">
        <td class="secretTableTop">Hammer Drop Test: {molecule.hammerDrop}</td>
        <td class="secretTableTop">Friction Test: {molecule.friction}</td>
        </tr>
        </table>
        <div class="pagebreak"> </div>
        <h3>Interpretation</h3>
        <p> {interpRoS} </p>
        <p> {interpOB}  </p>
        <p> {interpYoshida} </p>
        <p> {interpTd24} </p>
        <p> {interpOREOS} </p>
        <table align="Center" width="65%" >
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
        <td>{molecule.oreoSmallScale_val}</td>
        <td>{molecule.oreoTensScale_val}</td>
        <td>{molecule.oreoHundsScale_val}</td>
        <td>{molecule.oreoLargeScale_val}</td>
        </tr>
        </tbody>
        </table>
        <p><small>[1]: <i>Org. Proc. Res. Dev.,</i> 2021, 25, 2, 212-224</small><br>
        <small>[2]: <i>Org. Proc. Res. Dev.,</i> 2020, 24, 1, 67-84</small><br>
        <small>[3]: <i>Angew. Chem. Int. Ed.,</i> 2020, 59, 15798-15802</small></p>
        <div id="footer_content" style="text-align: center;">
        <p><small>This memo may contain confidential information.</small><br>
        <small><a href="https://github.com/Tomisall/ProcessSafteyDB">ThermalDex</a></small></p>
        </div>
        </body>'''
    
    if additionalData != '<h2>Additional Data</h2>':
        html_content += additionalData
    else:
        print('\n\nFoo\n\n')

    #now = datetime.now()
    #neatNow = now.strftime("%d-%b-%Y_%H-%M-%S")
    #outName = f'./_core/UserAddedData/temp/temp_{str(neatNow)}.pdf'
    convert_html_to_pdf(html_content, "./AssessmentMemos/ThermalHazardAssessmentMemo.pdf")
    time.sleep(1)

    if sys.platform != "darwin":
        subprocess.run(['ii', "./AssessmentMemos/ThermalHazardAssessmentMemo.pdf"], check=True, shell=True)

    else:
        subprocess.run(['open', "./AssessmentMemos/ThermalHazardAssessmentMemo.pdf"], check=True, shell=True)


def multiReportCreation(resultsTable):
    
    html_content = f'''
        <!DOCTYPE html>
        <html>
        <head>
        <META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=UTF-8">
        </head>
        <body>
        <p><link rel="stylesheet" href="./_core/style_multirep.css"></p>
        <div style="text-align: center;">
        <h1>Thermal Hazard Assessment Memo</h1>
        </div>
        {resultsTable}
        <p><small>[1]: <i>Org. Proc. Res. Dev.,</i> 2021, 25, 2, 212-224</small><br>
        <small>[2]: <i>Org. Proc. Res. Dev.,</i> 2020, 24, 1, 67-84</small><br>
        <small>[3]: <i>Angew. Chem. Int. Ed.,</i> 2020, 59, 15798-15802</small></p>
        <div id="footer_content" style="text-align: center;">
        <p><small>This memo may contain confidential information.</small><br>
        <small><a href="https://github.com/Tomisall/ProcessSafteyDB">ThermalDex</a></small></p>
        </div>
        </body>'''
    
    convert_html_to_pdf(html_content, "./AssessmentMemos/MultipleReport_ThermalHazardAssessmentMemo.pdf")
    time.sleep(1)

    if sys.platform != "darwin":
        subprocess.run(['start', "./AssessmentMemos/MultipleReport_ThermalHazardAssessmentMemo.pdf"], check=True, shell=True)

    else:
        subprocess.run(['open', "./AssessmentMemos/MultipleReport_ThermalHazardAssessmentMemo.pdf"], check=True, shell=True)