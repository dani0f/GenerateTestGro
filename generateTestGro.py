
class testGro:
    def __init__(self):
        self.dt = '0.1'
        self.population_max='2000'
        self.signals='1.0'
        self.signals_draw='1.0'
        self.iptg1=['1.0','0.0001'] #kdiff.kdeg
        self.ColagenoSig=['1.0','0.0001']
        self.ElpSig=['1.0','0.0001']
        self.ElpOperonActTimes='10.0,10.0,10.0'
        self.ElpOperonActVariabilities='4.5,6.6,6.6'
        self.ElpOperonDegTimes='10.0,10.0,10.0'
        self.ElpOperonDegVariabilities='1.5,2.0,2.0'
        self.ColagenoOperonActTimes='10.0,10.0,10.0'
        self.ColagenoOperonActVariabilities='4.5,6.6,6.6'
        self.ColagenoOperonDegTimes='6.0,10.0,10.0'
        self.ColagenoOperonDegVariabilities='1.5,2.0,2.0'
        self.UnionOperonActTimes='0.1'
        self.UnionOperonActVariabilities='0.0'
        self.UnionOperonDegTimes='0.0'
        self.UnionOperonDegVariabilities='0.0'
        self.iptg1Umbral='0.1' #> x ip1
        self.conjugateYFP='0.2' # p1
        self.conjugateGFP='0.2' # p2
        self.growthRate='0.017'
        self.ColagenoSigUmbral='10' #exact
        self.ElpSigUmbral='10' #exact
        self.ecolisP1='200,0,0,300'
        self.ecolisP2="200,0,0,300"

def generateCode(test): 
        code=""""  
        include gro
        set ( "dt", """+test.dt+""" );
        set ( "population_max", """+test.population_max+""" );
        set ("signals", """+test.signals +""");
        set ("signals_draw", """+test.signals_draw +""");
        grid("continuous", "gro_original", 10, 10, 8);
        
        iptg1 := s_signal([kdiff := """+test.iptg1[0]+""", kdeg:= """+test.iptg1[1]+"""]);
        ColagenoSig := s_signal([kdiff := """+test.ColagenoSig[0]+""", kdeg:= """+test.ColagenoSig[1]+"""]);
        ElpSig := s_signal([kdiff := """+test.ElpSig[0]+""", kdeg:= """+test.ElpSig[1]+"""]);

        
        genes([ name := "ElpOperon",
                proteins := {"Elp","YFP"},
                promoter := [function := "YES", transcription_factors := {"ip1"}],
                prot_act_times := [times := {"""+test.ElpOperonActTimes+"""}, variabilities := {"""+test.ElpOperonActVariabilities+"""}],
                prot_deg_times := [times := {"""+test.ElpOperonDegTimes+"""}, variabilities := {"""+test.ElpOperonDegVariabilities+"""}]
        ]);
        
        genes([ name := "ColagenoOperon",
                proteins := {"Colageno","GFP"},
                promoter := [function := "YES", transcription_factors := {"ip1"}],
                prot_act_times := [times := {"""+test.ColagenoOperonActTimes+"""}, variabilities := {"""+test.ColagenoOperonActVariabilities+"""}],
                prot_deg_times := [times := {"""+test.ColagenoOperonDegTimes+"""}, variabilities := {"""+test.ColagenoOperonDegVariabilities+"""}]
        ]);
        
        genes([ name := "UnionOperon",
                proteins := {"union"},
                promoter := [function := "AND", transcription_factors := {"GFP","YFP"}],
                prot_act_times := [times := {"""+ test.UnionOperonActTimes+"""}, variabilities := {"""+test.UnionOperonActVariabilities+"""}],
                prot_deg_times := [times := {"""+test.UnionOperonDegTimes+"""}, variabilities := {"""+test.UnionOperonDegVariabilities+"""}]
        ]);
        
        plasmids_genes([p1 := {"ElpOperon"},
                p2 := {"ColagenoOperon"},
                p3 := {"UnionOperon"}
        ]);
        //Actions
        action({"YFP"}, "d_paint", {"0","0","1","0"});
        action({"-YFP"}, "d_paint", {"0","0","-1","0"});
        action({"GFP"}, "d_paint", {"1","0","0","0"});
        action({"-GFP"}, "d_paint", {"-1","0","0","0"});
        action({"union"}, "paint", {"0","255","0","0"});
        action({"-union"}, "d_paint", {"0","-5","0","0"});
        action({}, "s_get_QS", {tostring(iptg1), ">", " """+test.iptg1Umbral +""" ", "ip1"});
        action({"YFP"}, "conjugate", {"p1"," """+test.conjugateYFP +""" "});
        action({"GFP"}, "conjugate", {"p2"," """+test.conjugateGFP +""" "});
        action({},"set_growth_rate",{" """+test.growthRate +""" "});
        action({"Colageno"},"s_emit_signal",{tostring(ColagenoSig)," """+test.ColagenoSigUmbral +""" ","exact"}); 
        action({"Elp"},"s_emit_signal",{tostring(ElpSig)," """+test.ElpSigUmbral +""" ","exact"}); 
        //Programs
        program p() :=
        {
                skip();
        };
        program main() := {
        true:
        {
                s_set_signal(iptg1, 10, 0, 0);
        }
        c_ecolis("""+test.ecolisP1 +""", {"p1","p3"}, program p());
        c_ecolis("""+test.ecolisP2 +""", {"p2","p3"}, program p());
        //c_ecolis(200, 0, 0, 300, {"p1","p2","p3"}, program p());
        
        };
        """        
        return(code)

test= testGro() #Objeto con parametros iniciales del código

#Ejemplo de creación de pruebas
#Genera 5 códigos con población maxima diferente
population=1000
for i in range(5):
        population=population+1000
        test.population_max = str(population)
        file = open("test"+str(i)+".gro", "w", encoding="utf-8")
        file.write(generateCode(test))
        file.close()