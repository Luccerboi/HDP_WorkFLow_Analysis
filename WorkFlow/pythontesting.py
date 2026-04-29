import os
from CsBBCl_ht_new import *

# os.chdir('/home/lucw/ht_programfiles/test')
# loc = initialize_compounds()
# p = ProcessLedger()
# #p.StartNewLedger(list_of_compounds = loc)


# #p.SetSingleValue(loc[2],"2a","completed",15)

# #p.SetSingleValue(loc[-1],4,"IncarExtras",{"EDIFF":"1E-8"})

# o = p.GetCompletionOverview()


# qd = p.GetQueue()

# print(o)

# print(qd)


def testWFcreation(TestJobs:str = "/home/lucw/ht_programfiles/TestJobs.json", CompList: str ="/home/lucw/ht_programfiles/testinput", TestDir: str = "/home/lucw/ht_programfiles/test", LedgerName:str = "JobInformationLedger.csv", cleanstart:bool=True):
    if os.path.exists(TestDir):
        os.chdir(TestDir)

        if cleanstart:
            os.system("rm -r *")
    else:
        os.mkdir(TestDir)
        os.chdir(TestDir)


    loc = initialize_compounds(CompList)
    ledger = ProcessLedger(JobSettingsPath=TestJobs, StartPath=TestDir,ledger_filename=LedgerName)
    if cleanstart:
        ledger.StartNewLedger(list_of_compounds=loc)

    print(ledger.GetCompletionOverview())


    LedgerDF = ledger.ReadFullLedger()
    print(LedgerDF.columns.to_list())

    print(ledger.GetQueue())

    return loc

def InsertPreviousPoscars(ListOfCompounds:list,ledger:ProcessLedger, NewDBDir:str, OldDBDir:str):

    for comp in ListOfCompounds:
        oldcompfolder = os.path.join(OldDBDir, repr(comp))
        if os.path.exists():
            pass


def UsePOTCARENMAX(MultFactor=1.2):

    #os.chdir(self.JobSpecifics["JobPath"])
    get_enmax = subprocess.run(["awk", "BEGIN{a = 0} /ENMAX/ {if($3>0+a) a=$3 fi} END{print a}", "POTCAR"], stdout=subprocess.PIPE, text=True)
    print(get_enmax.stdout)
    enout = float(get_enmax.stdout.strip().strip(';'))

    EMaxUse = MultFactor * enout

    return EMaxUse



if __name__ == "__main__":
    testWFcreation(cleanstart=True)
