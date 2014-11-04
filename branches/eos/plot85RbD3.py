#!/usr/bin/python

import ROOT as r
mycol = (r.kBlue, r.kRed, r.kGreen+2, r.kMagenta+1, r.kYellow-3, r.kCyan-6,
         r.kAzure, r.kOrange, r.kViolet, r.kSpring, r.kPink, r.kTeal-7)

nt = r.TNtuple("opData", "opData",
               "t:g0:g1:g2:g3:g4:f0:f1:f2:f3:f4:f5:f6:e0:e1:e2:e3:e4:e5:e6:e7:e8:e9:e10:e11:e12:e13:e14:e15:e16:e17:e18:e19:e20:e21:e22:e23:tot:pol:ali:exc");
nt.ReadFile("opData.dat")


gstate = []
gleg = r.TLegend(0.1, 0.7, 0.48, 0.9)
for i in range(5) :
    command = "g" + str(i) + ":t"
    n = nt.Draw(command, "", "goff")
    gstate.append(r.TGraph(int(n), nt.GetV2(), nt.GetV1()))
    gstate[i].GetYaxis().SetRangeUser(0, 1)
    gstate[i].SetLineColor(mycol[i])
    gstate[i].SetLineStyle((i+1)%10)
    gstate[i].SetLineWidth(3)
    gstate[i].SetTitle("|F=2,Mf=" + str(i-2))
    gleg.AddEntry(gstate[i], gstate[i].GetTitle(), "L")
    
gcan = r.TCanvas()
gstate[0].Draw("AL")
for i in range(1, 5) :
    gstate[i].Draw("LSAME")
gleg.Draw()

fstate = []
fleg = r.TLegend(0.1, 0.7, 0.48, 0.9)
for i in range(7) :
    command = "f" + str(i) + ":t"
    n = nt.Draw(command, "", "goff")
    fstate.append(r.TGraph(int(n), nt.GetV2(), nt.GetV1()))
    fstate[i].GetYaxis().SetRangeUser(0, 1)
    fstate[i].SetLineColor(mycol[i])
    fstate[i].SetLineStyle((i+1)%10)
    fstate[i].SetLineWidth(3)
    fstate[i].SetTitle("|F=3,Mf=" + str(i-3))
    fleg.AddEntry(fstate[i], fstate[i].GetTitle(), "L")

fcan = r.TCanvas()
fstate[0].Draw("AL")
for i in range(1, 7) :
    fstate[i].Draw("LSAME")
fleg.Draw()

estate = []
eleg = r.TLegend(0.1, 0.7, 0.48, 0.9)
for i in range(24) :
    command = "e" + str(i) + ":t"
    n = nt.Draw(command, "", "goff")
    estate.append(r.TGraph(int(n), nt.GetV2(), nt.GetV1()))

# F = 1
for i in range(0, 3) :
    estate[i].SetTitle("|F=1,Mf=" + str(i-1))

# F = 2
for i in range(3, 7) :
    estate[i].SetTitle("|F=2,Mf=" + str(i-5))

# F = 2
for i in range(7, 14) :
    estate[i].SetTitle("|F=3,Mf=" + str(i-10))

# F = 3
for i in range(14, 24) :
    estate[i].SetTitle("|F=4,Mf=" + str(i-18))

for i in range(len(estate)) :
    estate[i].GetYaxis().SetRangeUser(0, 1)
    estate[i].SetLineColor(mycol[i%len(mycol)])
    estate[i].SetLineStyle((i+1)%10)
    estate[i].SetLineWidth(3)
    eleg.AddEntry(estate[i], estate[i].GetTitle(), "L")
    

ecan = r.TCanvas()
estate[0].Draw("AL")
for i in range(1, len(estate)) :
    estate[i].Draw("LSAME")
eleg.Draw()
ecan.SetTitle("Excited States")

n = nt.Draw("exc:t", "", "goff")
photocan = r.TCanvas()
exc_total = r.TGraph(int(n), nt.GetV2(), nt.GetV1())
exc_total.SetTitle("Total ES Population")
exc_total.SetLineColor(mycol[0])
exc_total.SetLineWidth(3)
exc_total.Draw("AL")

print "Total excited state population at end: ",
print exc_total.GetY()[exc_total.GetN()-1]
raw_input("Enter to end")

