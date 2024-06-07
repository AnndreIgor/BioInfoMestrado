FROM ubuntu:latest

# Copia os arquivos
COPY requirements.txt /home/andre/BioPyhton

# Atualize os pacotes e instale o Python e o pip
RUN apt-get update && apt-get install -y \
    git \
    sudo \
    # python
    python3 \
    python3-pip \
    python3.10-venv \
    # BioPython
    clustalw \
    muscle \
    clustalo \
    mafft \
    probcons \
    t-coffee
    
# Adicione um novo usuário chamado "andre"
RUN useradd -m -s /bin/bash andre && \
    echo "andre:012345" | chpasswd && \
    usermod -aG sudo andre

# Exponha a porta do Jupyter Notebook
EXPOSE 8888

# Defina o diretório de trabalho
WORKDIR /home/andre

# Defina o usuário padrão para o contêiner
USER andre

# Cria um volume para a pasta /home/andre
VOLUME /home/andre

# Ativar venv e iniciar o Jupyter Notebook
CMD /bin/bash -c "source /home/andre/Biopython/venv/bin/activate \
                  && jupyter notebook --ip=0.0.0.0 --port=8888 --allow-root"


# docker build -t python-jupyter .
# docker run --name BioPyhton -it -p 8888:8888 -v Dados-Biopython:/home/andre python-jupyter

# docker start BioPyhton
# docker exec -it BioPyhton /bin/bash

# && pip install -r /home/andre/Biopython/requirements.txt \
                  


