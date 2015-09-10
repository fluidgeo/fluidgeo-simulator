#set -x

# Authors:
# Eduardo Garcia (bidu@lncc.br)
# Diego Volpatto (volpatto@lncc.br)

dirPadrao="dirExp00"
dirNew="dirExp01"

dirExp="$(pwd)"/${1:-${dirPadrao}}
dirExpNew="$(pwd)"/${2:-${dirNew}}

comandoRun="$(pwd)/rodarSimulador.sh $dirNew"
echo $comandoRun
eval $comandoRun

comando="diff $dirExp/disp.1 $dirExpNew/disp.1"
echo $comando
eval $comando

read -p "Correto? s ou n: " checkVar
 
#echo $checkVar;

if [ $checkVar = "s" ]; then 
	echo $checkVar
	fileName="$(date +%H_%M_%d_%m_%Yfontes)";
        echo criando o seguinte diretorio de copia $fileName;
	cp -r fontes $fileName;
	ls -ltr $fileName;
else
	exit 	
fi
