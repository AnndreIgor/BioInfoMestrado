import os
import time
import psutil
import threading
import pandas as pd
import shutil 

from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline, ClustalOmegaCommandline, MafftCommandline, ProbconsCommandline, TCoffeeCommandline

def alinhamento(INPUT_PATH,
                OUTPUT_PATH,
                INPUT_FILE, 
                algoritmo):

    OUTPUT_FILE = os.path.join(OUTPUT_PATH, INPUT_FILE.split('.')[0] + '_' + algoritmo + '.aln')

    # print(f' --- {algoritmo} ---')
    match algoritmo.lower():
        case 'clustalw':
            clustalw_cline = ClustalwCommandline("clustalw",
                                             infile=os.path.join(INPUT_PATH, INPUT_FILE),
                                             outfile=OUTPUT_FILE)
            stdout, stderr = clustalw_cline()
    
        case 'muscle':
            muscle_cline = MuscleCommandline(input=os.path.join(INPUT_PATH, INPUT_FILE),
                                             out=OUTPUT_FILE
                                             # ,clwout=True
                                            )
            
            stdout, stderr = muscle_cline()
            print(stdout)
    
        case 'clustalo':
            clustalo_cline = ClustalOmegaCommandline(
                infile=os.path.join(INPUT_PATH, INPUT_FILE),
                outfile=OUTPUT_FILE,
                verbose=True,
                auto=True
            )
            
            stdout, stderr = clustalo_cline()
    
        case 'mafft':
            mafft_cline = MafftCommandline(
                input=os.path.join(INPUT_PATH, INPUT_FILE)
            )
            
            stdout, stderr = mafft_cline()
            
            with open(OUTPUT_FILE, "w") as handle:
                handle.write(stdout)
    
        case 'probcons':
            probcons_cline = ProbconsCommandline(input=os.path.join(INPUT_PATH, INPUT_FILE),
                                                 clustalw=True)
            stdout, stderr = probcons_cline()
            
            with open(OUTPUT_FILE, "w") as handle:
                handle.write(stdout)
    
        case 't-coffee':
            tcoffee_cline = TCoffeeCommandline(infile=os.path.join(INPUT_PATH, INPUT_FILE),
                                           output="clustalw",
                                           outfile=OUTPUT_FILE)
            tcoffee_cline()

        case other:
            print('Algotitmo não implementado')

def construcao_arvore(INPUT_PATH,
                      OUTPUT_PATH,
                      OUTPUT_FORMAT,
                      sequencia_alinhada, 
                      modelo_evolutivo):

    nome = sequencia_alinhada.split('.')[0]
    arquivo_saida = f'{nome}_{modelo_evolutivo}.{OUTPUT_FORMAT}'

    try:
        with open(os.path.join(INPUT_PATH, sequencia_alinhada), "r") as file:
            alignment = AlignIO.read(file, "fasta")
    except ValueError:
        with open(os.path.join(INPUT_PATH, sequencia_alinhada), "r") as file:
            alignment = AlignIO.read(file, "clustal")
    
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor()
    
    # Calcula a matriz de distâncias
    dm = calculator.get_distance(alignment)
    
    match modelo_evolutivo.lower():
        case 'nj':
            tree = constructor.nj(dm)
        case 'upgma':
            tree = constructor.upgma(dm)
        case other:
            print('Modelo evolutivo não encontrado')

    Phylo.write(tree,
                os.path.join(
                    OUTPUT_PATH,
                    arquivo_saida),
                    OUTPUT_FORMAT)

def monitor_system():
    initial_disk_io = psutil.disk_io_counters()
    while monitorar:
        # Utilização de CPU
        cpu_usage = psutil.cpu_percent(interval=1)
        memory_info = psutil.virtual_memory()
        disk_usage = psutil.disk_usage('/')
        net_io = psutil.net_io_counters()
        current_disk_io = psutil.disk_io_counters()
        read_bytes = current_disk_io.read_bytes - initial_disk_io.read_bytes
        write_bytes = current_disk_io.write_bytes - initial_disk_io.write_bytes

        dados_monitoramento = {
            'Tarefa': tarefa,
            'Algoritmo': algoritmo,
            'Arquivo de Entrada': arquivo_entrada,
            'Modelo Evolutivo': modelo_evolutivo,
            'Formato Saida': formato_saida,
            'Utilização de CPU': cpu_usage,
            'Memória total:': memory_info.total / (1024 ** 3),
            'Memória disponível': memory_info.available / (1024 ** 3),
            'Uso de memória': memory_info.percent,
            'Espaço total em disco': disk_usage.total / (1024 ** 3),
            'Espaço usado em disco': disk_usage.used / (1024 ** 3),
            'Espaço livre em disco': disk_usage.free / (1024 ** 3),
            'Uso de disco': disk_usage.percent,
            'Bytes lidos': read_bytes / (1024 ** 2),
            'Bytes escritos': write_bytes / (1024 ** 2),
            'Bytes enviados': net_io.bytes_sent / (1024 ** 2),
            'Bytes recebidos': net_io.bytes_recv / (1024 ** 2),
            'Timestamp': time.time()
        }
        dados.append(dados_monitoramento)
        time.sleep(5)  # Intervalo de 5 segundos

def limpar_saidas():
    shutil.rmtree(os.path.join('files', 'output'))
    os.mkdir(os.path.join('files', 'output'))
    os.mkdir(os.path.join('files', 'output', 'arvores_filogeneticas'))
    os.mkdir(os.path.join('files', 'output', 'sequencias_alinhadas'))


limpar_saidas()

monitorar = True
dados = []

tarefa = None
arquivo_entrada = None
modelo_evolutivo = None
formato_saida = None
algoritmo = None

# Inicia o monitor do sistema
monitor_thread = threading.Thread(target=monitor_system)
monitor_thread.daemon = True  # Permite que a thread seja finalizada quando o programa principal termina
monitor_thread.start()

INPUT_PATH = os.path.join('files', 'input')
OUTPUT_PATH = os.path.join('files', 'output', 'sequencias_alinhadas')
INPUT_FILE = 'ls_orchid.fasta'

algoritmos = ['clustalw', 'muscle', 'clustalo', 'mafft', 'probcons', 't-coffee']

tarefa = 'Sequenciamento'
for algoritmo in algoritmos:
    arquivo_entrada = INPUT_FILE
    alinhamento(INPUT_PATH,
                OUTPUT_PATH,
                INPUT_FILE,
                algoritmo)

INPUT_PATH = os.path.join('files', 'output', 'sequencias_alinhadas')
OUTPUT_PATH = os.path.join('files', 'output', 'arvores_filogeneticas')

modelos_evolutivos = ['nj', 'upgma']
formatos_saida = ['nexus', 'newick']

tarefa = 'Criação de Árvore'
for sequencia_alinhada in os.listdir(INPUT_PATH):
    if sequencia_alinhada.endswith('.aln'):
        for modelo_evolutivo in modelos_evolutivos:
            for formato_saida in formatos_saida:
                if sequencia_alinhada.endswith('.aln'):
                    arquivo_entrada = sequencia_alinhada
                    construcao_arvore(INPUT_PATH,
                                      OUTPUT_PATH,
                                      formato_saida,
                                      sequencia_alinhada,
                                      modelo_evolutivo)

# Para o monitoramento
monitorar = False
aux = pd.DataFrame(dados)
aux.to_csv('saida.csv', index=False) {:.2f}