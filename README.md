# Projeto Metodos Numericos:
Implementação dos metodos numericos: Euler,Euler Inverso, Euler Aprimorado, Runge Kutta, Adam-Bashforth, Adam-Multon e Formula Inversa.

# Pré-requisitos:
1.Instale a biblioteca sympy com os comandos:
  sudo apt-get install python3-pip
  sudo pip3 install sympy
  
2. Instale a biblioteca matplotlib:
  sudo apt-get install python3-matplotlib

3. Já no diretorio onde está o arquivo projeto.py execute:
  python3 projeto.py
 
# Executando o projeto.py:
 Ao executar o arquivo, será feito a analise de todos os exemplos de testes disponibilizados na especificação do projeto além dos que foram inseridos para ser analisados outros quesitos, como os metodoso calculados sem previsão, a cada execução porém, será possivel visualizar o gráfico da função, então, para que o proximo método seja executado é preciso fechar a janela do gráfico.
 
 
# Testes:
 Os testes utilizados na analise dos métodos estão no arquivo.txt, os valores calculados em cada método serão salvos no arquivo saida.txt, de acordo com a especificação do projeto, estes são alguns exemplos do tipo e formato de informação que será encontrado no arquivo.txt:
 
* euler 0 0 20 1-t+4*y
* euler_inverso 0 0 0.1 20 1-t+4*y
* euler_aprimorado 0 0 0.1 20 1-t+4*y
-runge_kutta 0 0 0.1 20 1-t+4*y
-adam_bashforth 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
-adam_multon 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
-formula_inversa 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6

Para os casos em que irá se analisar os metodos (bônus) sem previsão, foi inserido o nome _sem_previsao, na frente do nome de seu metodo, exemplo:
euler_inverso_sem_previsao 0 0 0.1 20 1-t+4*y
euler__aprimorado_sem_previsao 0 0 0.1 20 1-t+4*y

Também foi relizado a implementação (bônus) do método de Runge-Kutta de 5 ordem, que estará informado no aquivo.txt da seguinte forma:
runge_kutta5 0 0 0.1 20 1-t+4*y



 
