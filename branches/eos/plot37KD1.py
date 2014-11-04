#!/usr/bin/python

import ROOT as r
mycol = (r.kBlue, r.kRed, r.kGreen+2, r.kMagenta+1, r.kYellow-3, r.kCyan-6,
         r.kAzure, r.kOrange, r.kViolet, r.kSpring, r.kPink, r.kTeal-7)

nt = r.TNtuple("opData", "opData",
               "t:g0:g1:g2:f0:f1:f2:f3:f4:e0:e1:e2:e3:e4:e5:e6:e7:tot:pol:ali:exc");
nt.ReadFile("opData.dat")


gstate = []
gleg = r.TLegend(0.1, 0.7, 0.48, 0.9)
for i in range(3) :
    command = "g" + str(i) + ":t"
    n = nt.Draw(command, "", "goff")
    gstate.append(r.TGraph(int(n), nt.GetV2(), nt.GetV1()))
    gstate[i].GetYaxis().SetRangeUser(0, 1)
    gstate[i].SetLineColor(mycol[i])
    gstate[i].SetLineStyle((i+1)%10)
    gstate[i].SetLineWidth(3)
    gstate[i].SetTitle("|F=1,Mf=" + str(i-1))
    gleg.AddEntry(gstate[i], gstate[i].GetTitle(), "L")
    
gcan = r.TCanvas()
gstate[0].Draw("AL")
gstate[1].Draw("LSAME")
gstate[2].Draw("LSAME")
gleg.Draw()

fstate = []
fleg = r.TLegend(0.1, 0.7, 0.48, 0.9)
for i in range(5) :
    command = "f" + str(i) + ":t"
    n = nt.Draw(command, "", "goff")
    fstate.append(r.TGraph(int(n), nt.GetV2(), nt.GetV1()))
    fstate[i].GetYaxis().SetRangeUser(0, 1)
    fstate[i].SetLineColor(mycol[i])
    fstate[i].SetLineStyle((i+1)%10)
    fstate[i].SetLineWidth(3)
    fstate[i].SetTitle("|F=2,Mf=" + str(i-2))
    fleg.AddEntry(fstate[i], fstate[i].GetTitle(), "L")

fcan = r.TCanvas()
fstate[0].Draw("AL")
for i in range(1, 5) :
    fstate[i].Draw("LSAME")
fleg.Draw()

print "Minimum |2,2> population ",
print fstate[4].GetY()[r.TMath.LocMin(fstate[4].GetN(), fstate[4].GetY())]

estate = []
eleg = r.TLegend(0.1, 0.7, 0.48, 0.9)
for i in range(8) :
    command = "e" + str(i) + ":t"
    n = nt.Draw(command, "", "goff")
    estate.append(r.TGraph(int(n), nt.GetV2(), nt.GetV1()))

# F = 1
for i in range(0, 3) :
    estate[i].SetTitle("|F=1,Mf=" + str(i-2))

# F = 2
for i in range(3, 8) :
    estate[i].SetTitle("|F=2,Mf=" + str(i-6))

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
exc_total.SetNameTitle("ES Population", "ES Population")
exc_total.SetLineColor(mycol[0])
exc_total.SetLineWidth(3)
exc_total.Draw("AL")
print "Final E.S. population ", exc_total.GetY()[exc_total.GetN()-1]



n = nt.Draw("pol:t", "", "goff")
pcan = r.TCanvas()
pol = r.TGraph(int(n), nt.GetV2(), nt.GetV1())
pol.SetNameTitle("Polarization", "Polarization")
pol.SetLineColor(mycol[0])
pol.SetLineWidth(3)
pol.Draw("AL")

print "Final polarization: ", pol.GetY()[pol.GetN()-1]

raw_input("Enter to end")

