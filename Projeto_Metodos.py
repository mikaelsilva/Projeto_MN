from sympy import symbols, Function, solve,Eq,S
import matplotlib.pyplot as plt
import sympy 




def euler(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y,lista_t):
	
	print ('Metodo de euler')
	print ('y(',valor_t,') =',valor_y)
	print ('h =',valor_h)
	lista_y.append(valor_y)
	lista_t.append(valor_t)

	for next in range(valor_passos+1):
		Fn=S(expressao).subs({'t':valor_t,'y':valor_y})
		valor_y=valor_y+valor_h*Fn
		valor_t=valor_t+valor_h

		lista_y.append(valor_y)
		lista_t.append(valor_t)	
		print ("%d  %.10f" %(next , valor_y))

	
	plt.plot(lista_t, lista_y)
	plt.show()
	return 

def euler_n_mais_1(valor_y2,valor_t2,valor_h2):

	Fn2=S(expressao).subs({'t':valor_t2,'y':valor_y2})
	valor_y2=valor_y2+valor_h2*Fn2
	
	return valor_y2

def euler_inverso(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_inverso,lista_t_inverso):

	print ('Metodo de euler_inverso')
	print ('y(',valor_t,') =',valor_y)
	print ('h =',valor_h)
	lista_y_inverso.append(valor_y)
	lista_t_inverso.append(valor_t)
	
	for next in range(valor_passos):
		
		valor_y_mais_1=euler_n_mais_1(valor_y,valor_t,valor_h) #Yn+1 para Fn+1
		valor_t_mais_1=valor_t+valor_h #Tn+1

		Fn=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1}) #Fn+1
		valor_y=valor_y+valor_h*Fn
		valor_t=valor_t+valor_h
		print ("%.1f  %.10f" %(valor_t , valor_y))
		lista_y_inverso.append(valor_y)
		lista_t_inverso.append(valor_t)
	return 

def euler_aprimorado(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_aprimorado,lista_t_aprimorado):

	print ('Metodo de euler_aprimorado')
	print ('y(',valor_t,') =',valor_y)
	print ('h =',valor_h)
	lista_y_aprimorado.append(valor_y)
	lista_t_aprimorado.append(valor_t)
	
	for next in range(valor_passos):

		#Fn_mais_1
		valor_y_mais_1=euler_n_mais_1(valor_y,valor_t,valor_h)
		valor_t_mais_1=valor_t+valor_h
		Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

		#Fn
		Fn=S(expressao).subs({'t':valor_t,'y':valor_y})
		Fn_total=Fn+Fn_mais_1

		valor_y=valor_y+(valor_h/2)*Fn_total
		valor_t=valor_t+valor_h

		print ("%.1f  %f" %(valor_t , valor_y))
		lista_y_aprimorado.append(valor_y)
		lista_t_aprimorado.append(valor_t)

	return

def rugge_kutta(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_kutta,lista_t_kutta):

	print ('Metodo de rugge_kutta')
	print ('y(',valor_t,') =',valor_y)
	print ('h =',valor_h)
	lista_y_kutta.append(valor_y)
	lista_t_kutta.append(valor_t)
	
	for next in range(valor_passos):

		K_1=S(expressao).subs({'t':valor_t,'y':valor_y})
		K_2=S(expressao).subs({'t': valor_t+(0.5)*valor_h ,'y':valor_y+(0.5*valor_h*K_1)})
		K_3=S(expressao).subs({'t': valor_t+(0.5)*valor_h ,'y':valor_y+(0.5*valor_h*K_2)})
		K_4=S(expressao).subs({'t': valor_t+valor_h ,'y':valor_y+valor_h*K_3})

		Fn=(K_1+ 2*K_2 + 2*K_3 + K_4)
		valor_y=valor_y+(valor_h/6)*Fn
		valor_t=valor_t+valor_h

		print ("%.1f  %f" %(valor_t , valor_y))
		lista_y_kutta.append(valor_y)
		lista_t_kutta.append(valor_t)

	return 

def rugge_kutta5(expressao,valor_y,valor_t,valor_h,valor_passos,valor_y_kutta5,valor_t_kutta5):
	
	print ('Metodo de rugge_kutta 5 Ordem')
	print ('y(',valor_t,') =',valor_y)
	print ('h =',valor_h)
	lista_y_kutta5.append(valor_y)
	lista_t_kutta5.append(valor_t)
	
	for next in range(valor_passos):

		K_1=S(expressao).subs({'t':valor_t,'y':valor_y})
		
		T_k_2=valor_h/4
		Y_k_2=((valor_h/4)*valor_h*K_1)
		K_2=S(expressao).subs({'t': valor_t+T_k_2,'y':valor_y+Y_k_2})

		T_k_3=valor_h/4
		Y_k_3=(valor_h*K_1)/8+(K_2*valor_h)/8
		K_3=S(expressao).subs({'t': valor_t+T_k_3,'y':valor_y+Y_k_3})

		T_k_4=valor_h/2
		Y_k_4=(K_3*valor_h)-(valor_h*K_2)/2
		K_4=S(expressao).subs({'t': valor_t+T_k_4,'y':valor_y+Y_k_4})

		T_k_5=(valor_h*(3/4))
		Y_k_5=(K_1*valor_h)*(3/16)+((valor_h*K_4)*(9/16))
		K_5=S(expressao).subs({'t': valor_t+T_k_5,'y':valor_y+Y_k_5})

		Y_k_6=(K_2*valor_h)*(2/7)+((valor_h*K_3)*(12/7))-((3/7)*K_1*valor_h)-((12/7)*K_4*valor_h)+((8/7)*K_5*valor_h)
		K_6=S(expressao).subs({'t': valor_t+valor_h,'y':valor_y+Y_k_5})


		Fn=(7*K_1 + 32*K_3 + 12*K_4 + 7*K_6)

		valor_y=valor_y+(valor_h/90)*Fn
		valor_t=valor_t+valor_h

		print ("%.1f  %f" %(valor_t , valor_y))
		lista_y_kutta.append(valor_y)
		lista_t_kutta.append(valor_t)
	return 

def adam_bashforth(expressao,leitura,lista_y_adam,lista_t_adam):
	ordem='2'
	if(ordem=='2'):
		valor_h=float(leitura[4])
		valor_passos=int(leitura[5])
		valor_t_inicial=lista_t_adam[0]
		valor_y=lista_y_adam[1]
		valor_y_menos_1=lista_y_adam[0]
		valor_t=2*valor_h

		print (valor_t_inicial,valor_y_menos_1)
		print (valor_t,valor_y)

		for next in range(valor_passos):	
			Fn_menos_1=S(expressao).subs({'t':valor_t_inicial,'y':valor_y_menos_1})
			Fn=S(expressao).subs({'t':valor_t,'y':valor_y})
			valor_y=valor_y+((3/2)*valor_h*Fn)-((0.5)*valor_h*Fn_menos_1)
			
			valor_t_inicial=valor_t_inicial+valor_h
			valor_t=valor_t+valor_h

			print ("%.1f %f" %(valor_t , valor_y))

			lista_t_adam.append(valor_t)
			lista_y_adam.append(valor_y)

		print ('Terminado')

	return 


#Lendo o arquivo
arquivo=open('arquivo.txt','r')
leitura=arquivo.read().split()


#Declarações Iniciais
valor_x=0
valor_y=0
valor_h=0
valor_passos=0



#Guardando os valores
#valor_t=float(leitura[1])
#valor_y=float(leitura[2])
#valor_h=float(leitura[3])
#valor_passos=int(leitura[4])
#expressao=leitura[5]


#Listas Usadas
lista_y=[]
lista_t=[]

lista_y_inverso=[]
lista_t_inverso=[]

lista_y_aprimorado=[]
lista_t_aprimorado=[]

lista_y_kutta=[]
lista_t_kutta=[]

lista_y_kutta5=[]
lista_t_kutta5=[]

lista_y_adam=[]
lista_t_adam=[]

#euler(valor_y,valor_t,valor_h,valor_passos,lista_y,lista_t)
#euler_inverso(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_inverso,lista_t_inverso)
#euler_aprimorado(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_aprimorado,lista_t_aprimorado)
#rugge_kutta(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_kutta,lista_t_kutta)
#rugge_kutta5(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_kutta5,lista_t_kutta5)
#plt.plot(lista_t_kutta, lista_y_kutta)
#plt.show()


if(leitura[0]=='adam_bashforth'):
	if  (leitura[len(leitura)-1]=='2'):
		y1=float(leitura[1])
		lista_y_adam.append(y1)
		y1=float(leitura[2])
		lista_y_adam.append(y1)
		y1=float(leitura[3])
		lista_t_adam.append(y1)
		adam_bashforth(leitura[6],leitura,lista_y_adam,lista_t_adam)
	elif(leitura[len(leitura)-1]=='3'):
	elif(leitura[len(leitura)-1]=='4'):
	elif(leitura[len(leitura)-1]=='5'):
	elif(leitura[len(leitura)-1]=='6'):
	elif(leitura[len(leitura)-1]=='7'):
	elif(leitura[len(leitura)-1]=='8'):




#plt.plot(lista_t_adam, lista_y_adam)
#plt.show()


arquivo.close()

