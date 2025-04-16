import sys
import subprocess
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table #, TableStyles
from thermDex.thermDexMolecule import thermalDexMolecule
from dataclasses import asdict
from textwrap import wrap

#arial_path = "/System/Library/Fonts/Supplemental/Arial.ttf"
#arial = TTFont("Arial","/System/Library/Fonts/Supplemental/Arial.ttf")
pdfmetrics.registerFont(TTFont('Arial', '.\\_core\\Arial.ttf'))
memoFont = "Arial"

pdfmetrics.registerFont(TTFont('Consolas', '.\\_core\\consola.ttf'))
unicodeFont = "Consolas"

my_Style=ParagraphStyle('My Para style',
fontName='Times-Roman',
backColor='#F1F1F1',
fontSize=16,
borderColor='#FFFF00',
borderWidth=2,
borderPadding=(20,20,20),
leading=20,
alignment=0
)

def apply_scripting(textobject, text, rise):
    textobject.setFont(memoFont, 8)
    textobject.setRise(rise)
    textobject.textOut(text)
    textobject.setFont(memoFont, 11)
    textobject.setRise(0)


def create_pdf(Name, results, image_path):
    filename = "./AssessmentMemos/thermal_assessment_memo.pdf"

    # Create the PDF
    c = canvas.Canvas(filename, pagesize=A4)
    c.setAuthor("ThermalDex - A thermal hazard assessment software")

    fontsAva = c.getAvailableFonts()
    print(fontsAva)

    # Add title
    c.setFont(memoFont, 16)
    c.drawCentredString(300, 775, "Thermal Hazard Assessment Memo")
    
    # Name Molecule
    c.setFont(memoFont, 14)
    c.drawString(50, 725, Name)

    # Add image at the top center
    if image_path:
        c.drawInlineImage(image_path, 147.5, 555, width=300, height=150) #for 297.5 centre A4 of 225 width img pos = 185, 500

    startOfText = 515
    # Add properties section
    c.setFont(memoFont, 11)
    c.drawString(50, startOfText, "Properties:")

    # Add the actual properties from your assessment
    wraped_SMILES = "\n".join(wrap(results["SMILES"], 65))
    textobject = c.beginText()
    textobject.setTextOrigin(70, startOfText-20)
    textobject.textLines("SMILES: " + wraped_SMILES)   
    c.drawText(textobject)
    #c.writeLines(70, startOfText-20, "SMILES: " + wraped_SMILES)
    try:
        c.drawString(70, startOfText-50, "Name: " + results["name"])
    except:
        c.drawString(70, startOfText-50, "Name: ")
    try:
        c.drawString(70, startOfText-70, "Formula: " + results["eleComp"])
    except:
        c.drawString(70, startOfText-70, "Formula: ")   
    try:
        c.drawString(70, startOfText-90, "mp: " + str(results["mp"]) + " to " +  str(results["mpEnd"] + " °C"))
#        textobject = c.beginText()
#        textobject.setTextOrigin(200, 375)
#        textobject.textOut('37')
#        apply_scripting(textobject, '%', -4)
#        c.drawText(textobject)
    except:
        c.drawString(70, startOfText-90, "mp: " + " " + " to " + " ")

    # Add results section
    c.setFont(memoFont, 11)
    c.drawString(50, startOfText-120, "Results:")

    # Add the actual results from your assessment
    # for i, (key, value) in enumerate(results.items()):
    #    c.drawString(70, 315 - i * 20, f"{key}: {value}")
    c.drawString(70, startOfText-140, "High Energy Groups =  " + str(results["HEG"]) + " "  + ", ".join(results["HEG_list"]))
    c.drawString(70, startOfText-160, "Explosive  Groups =  " + str(results["EFG"]) + " " + ", ".join(results["EFG_list"]))

    textobject = c.beginText()
    textobject.setTextOrigin(70, startOfText-200)
    textobject.textOut("Q")   
    apply_scripting(textobject, "DSC", -4)
    
    try:
        textobject.textOut(" = " + "{:.2f}".format(results["Q_dsc"]))
        textobject.setFont(unicodeFont, 11)
        textobject.textOut(" " + str(results["Qunits"]))
        textobject.setFont(memoFont, 11)
        c.drawText(textobject)
    except:
        c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(230, startOfText-200)
    textobject.textOut("T")   
    apply_scripting(textobject, "onset", -4)

    try:
        textobject.textOut(" = " + "{:.2f}".format(results["onsetT"]))
        textobject.textOut(" °C")
        c.drawText(textobject)
    except:
        c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(430, startOfText-200)
    textobject.textOut("T")   
    apply_scripting(textobject, "init", -4)
  
    try:
        textobject.textOut(" = " + "{:.2f}".format(results["initT"]))
        textobject.textOut(" °C")
        c.drawText(textobject)
    except:
        c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(70, startOfText-180)
    textobject.textOut("Rule of Six")   
    #apply_scripting(textobject, "init", -4)
    textobject.textOut(" = " + str(results["RoS_val"]))
    c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(230, startOfText-180)
    textobject.textOut("Oxygen Balance")   
    #apply_scripting(textobject, "init", -4)
    textobject.textOut(" = " + "{:.2f}".format(results["OB_val"]))
    c.drawText(textobject)


    textobject = c.beginText()
    textobject.setTextOrigin(70, startOfText-220)
    textobject.textOut("Impact Sensitivity")   

    try:
        #apply_scripting(textobject, "init", -4)
        textobject.textOut(" = " + "{:.2f}".format(results["IS_val"]))
        c.drawText(textobject)
    except:
        c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(230, startOfText-220)
    textobject.textOut("Explosive Propagation")   
    
    try:
        #apply_scripting(textobject, "init", -4)
        textobject.textOut(" = " + "{:.2f}".format(results["EP_val"]))
        c.drawText(textobject)
    except:
        c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(430, startOfText-220)
    textobject.textOut("T")   
    apply_scripting(textobject, "D24", -4)
    
    try:
        #d24Str = "{:.1f}".format(d24Temp)
        textobject.textOut(" = " + "{:.1f}".format(results["Td24"]))
        textobject.textOut(" °C")
        c.drawText(textobject)
    except:
        c.drawText(textobject)

    oreoSmall = results["oreoSmallScale_des"]
    oreoTens = results["oreoTensScale_des"]
    oreosHund = results["oreoHundsScale_des"]
    oreosLarge = results["oreoLargeScale_des"]
    oreosData = [["<5 g", "5 to 100 g", "100 to 500 g", ">500 g"], [oreoSmall, oreoTens, oreosHund, oreosLarge]]
    oreosTable = Table(oreosData, style=[
                ('GRID',(0,0),(-1,-1),0.5,colors.grey),
                #('BACKGROUND',(0,0),(1,1),colors.palegreen),
                #('SPAN',(0,0),(1,1)),
                #('BACKGROUND',(-2,-2),(-1,-1), colors.pink),
                #('SPAN',(-2,-2),(-1,-1)),
                ])
    oreosTable.wrapOn(c, 0, 0)
    oreosTable.drawOn(c, 160, startOfText-290)

    # Add interpretation section
    c.drawString(50, startOfText-325, "Interpretation:")
    c.setFont(memoFont, 11)
    # Add your interpretation text
    #interpretation_text = interpretation
    # interpPara = Paragraph(interpretation_text, my_Style)
    # interpPara.drawOn(c, 70, 480)

    # i = startOfText-345
    # for line in interpretation_text:
    #    c.drawString(70, i, line)
    #    i -= 20

    textobject = c.beginText()
    textobject.setTextOrigin(70, startOfText-345)
    #wraped_text = "\n".join(wrap(text, 80))
    textobject.textLines("The Rule of Six¹ value imples" + str(results["RoS_des"]) + ". The Oxygen Balance¹ suggests " + str(results["OB_des"]) + ".\nThe " + str(results["yoshidaMethod"]) + " method² was used to calculate Impact Sensitivity and Explosive Propagation values, these\n suggest " + str(results["IS_des"]) + " and " + str(results["EP_des"]) + ".")   
    c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(70, startOfText-385)
    textobject.textOut("The T")   
    apply_scripting(textobject, "D24", -4)
    #d24Str = "{:.1f}".format(d24Temp)
    textobject.textOut(" result gives the maximum safe operation temperature.")
    #textobject.textOut(" °C")
    c.drawText(textobject)


    # Add footer with disclaimer
    c.setFont(memoFont, 8)
    c.drawString(50, 70, "[1] Org. Proc. Res. Dev., 2011, 2341-2356")
    c.drawString(50, 60, "[2] Org. Proc. Res. Dev., 2021, 2117-2119")
    disclaimer_text = "This report may contain confidential information."
    c.drawCentredString(300, 30, disclaimer_text)

    # Save the PDF
    c.save()
    print(f"PDF report generated: {filename}")
    # `cwd`: current directory is straightforward
    #cwd = Path.cwd()
    cwd = "./"
    filepath =  cwd + filename #(cwd / filename).resolve()
    print(filepath)
    subprocess.run(['start', filepath], check=True, shell=True)

# Example usage:
#molecule = thermalDexMolecule(SMILES='CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',  name='Sildenafil very very very very very very #long name', Q_dsc=770.26, Qunits='J g⁻¹', onsetT=90.14, initT=74.62, proj='PDF_Test') 
#molecule = thermalDexMolecule(SMILES='CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', name='TNT', Q_dsc=770.26, Qunits='J g⁻¹', onsetT=90.14, #initT=74.62, proj='PDF_Test')
#molecule.genAllValues()
#Name = molecule.name
#results = asdict(molecule) 
#image_path = molecule.molIMG #"./_core/ThermalDexIcon.jpg"

#create_pdf(Name, results, image_path)
