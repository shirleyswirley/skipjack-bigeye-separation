# - Month fractions corresp to jan, feb, mar, etc.
monthfracs = [0,0.0834,0.1666,0.25,0.3334,0.4166,0.5,0.5834,0.6666,0.75,0.8334,0.9166];

# - Load ONI 
df = pd.read_csv('/ltraid2/sleung/data/ENSOindices/oni/oniindex1950_2018.txt',names=['Date','ONI']);
#print(df)
onifulltime = df['Date'];
onibegmoidx = np.where(onifulltime == begdate.year[0] + monthfracs[begdate.month[0]-1])[0][0];
oniendmoidx = np.where(onifulltime == enddate.year[0] + monthfracs[enddate.month[0]-1])[0][0];
onifull = df['ONI'];
oni = onifull[onibegmoidx:oniendmoidx+1].as_matrix();
onien = [0]*len(oni) # 1 = el nino month
oniln = [0]*len(oni) # 1 = la nina month

cmcounter = 0; # consecutive months counter
for imonth in range(len(oni)):
    if oni[imonth]>=0.5:
        cmcounter=cmcounter+1;
    elif oni[imonth]<0.5:
        cmcounter=0;
    if cmcounter>=5:
        onien[imonth-cmcounter:imonth]=[1]*cmcounter;

cmcounter = 0; # consecutive months counter
for imonth in range(len(oni)):
    if oni[imonth]<=-0.5:
        cmcounter=cmcounter+1;
    elif oni[imonth]>-0.5:
        cmcounter=0;
    if cmcounter>=5:
        oniln[imonth-cmcounter:imonth]=[1]*cmcounter;
        
if plotensoidx==1:
    fig = plt.figure(figsize=(11,4));
    plt.plot(onitime,oni);
    plt.plot(onitime,onien);
    plt.plot(onitime,oniln);
    plt.legend(["oni","onien","oniln"]);
    fig.show();
