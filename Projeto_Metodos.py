from sympy import symbols, Function, solve,Eq,S
import matplotlib.pyplot as plt
import sympy 



def euler(expressao,valor_h,valor_passos,lista_y_euler,lista_t_euler,euler):
	
	valor_y=lista_y_euler[0]
	valor_t=lista_t_euler[0]
	if(euler==1):print (valor_t,  valor_y)

	for i in range(valor_passos):
		Fn=S(expressao).subs({'t':lista_t_euler[i],'y':lista_y_euler[i]})
		valor_y=valor_y+valor_h*Fn
		valor_t=valor_t+valor_h

		lista_y_euler.append(valor_y)
		lista_t_euler.append(valor_t)	
		if(euler==1):print ("%d  %.10f" %(i+1 , valor_y))

	return 

def euler_n_mais_1(expressao,valor_y2,valor_t2,valor_h2):

	Fn2=S(expressao).subs({'t':valor_t2,'y':valor_y2})
	valor_y2=valor_y2+valor_h2*Fn2
	
	return valor_y2

def euler_inverso(expressao,valor_h,valor_passos,lista_y_inverso,lista_t_inverso,inverso):
	
	valor_y=lista_y_inverso[0]
	valor_t=lista_t_inverso[0]
	if(inverso==1):print (valor_t,  valor_y)

	
	for i in range(valor_passos):
		
		valor_y_mais_1=euler_n_mais_1(expressao,valor_y,valor_t,valor_h) #Yn+1 para Fn+1
		valor_t_mais_1=valor_t+valor_h #Tn+1

		Fn=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1}) #Fn+1
		valor_y=valor_y+valor_h*Fn
		valor_t=valor_t+valor_h
		if(inverso==1):print ("%d  %.10f" %(i , valor_y))
		
		lista_y_inverso.append(valor_y)
		lista_t_inverso.append(valor_t)
	return 

def euler_aprimorado(expressao,valor_h,valor_passos,lista_y_aprimorado,lista_t_aprimorado,aprimorado):

	valor_y=lista_y_aprimorado[0]
	valor_t=lista_t_aprimorado[0]
	if(aprimorado==1):print (valor_t,  valor_y)
	
	for i in range(valor_passos):

		#Fn_mais_1
		valor_y_mais_1=euler_n_mais_1(expressao,valor_y,valor_t,valor_h)
		valor_t_mais_1=valor_t+valor_h
		Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

		#Fn
		Fn=S(expressao).subs({'t':lista_t_aprimorado[i],'y':lista_y_aprimorado[i]})
		Fn_total=Fn+Fn_mais_1

		valor_y=valor_y+(valor_h/2)*Fn_total
		valor_t=valor_t+valor_h

		if (aprimorado==1): print ("%d  %f" %(i, valor_y))
		lista_y_aprimorado.append(valor_y)
		lista_t_aprimorado.append(valor_t)

	return

def runge_kutta(expressao,valor_h,valor_passos,lista_y_kutta,lista_t_kutta,kutta):


	valor_y=lista_y_kutta[0]
	valor_t=lista_t_kutta[0]
	if(kutta==1): print (valor_t,  valor_y)

	for i in range(valor_passos):

		K_1=S(expressao).subs({'t':valor_t,'y':valor_y})
		K_2=S(expressao).subs({'t': valor_t+(0.5)*valor_h ,'y':valor_y+(0.5*valor_h*K_1)})
		K_3=S(expressao).subs({'t': valor_t+(0.5)*valor_h ,'y':valor_y+(0.5*valor_h*K_2)})
		K_4=S(expressao).subs({'t': valor_t+valor_h ,'y':valor_y+valor_h*K_3})

		Fn=(K_1+ 2*K_2 + 2*K_3 + K_4)
		valor_y=valor_y+(valor_h/6)*Fn
		valor_t=valor_t+valor_h

		if(kutta==1): print ("%d  %f" %(i, valor_y))
		lista_y_kutta.append(valor_y)
		lista_t_kutta.append(valor_t)

	return 

def rumge_kutta5(expressao,valor_y,valor_t,valor_h,valor_passos,valor_y_kutta5,valor_t_kutta5):
	
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




#BASHFORTH 2,3,4,5.....
def adam_bashforth2(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam):
	Fn=[]
	for j in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_adam[n],'y':lista_y_adam[n]}))
		n=n+1
	Fn_0=(3*Fn[1])
	Fn_1=(1*Fn[0])
	Fn_total=((valor_h/2)*(Fn_0-Fn_1))
	return Fn_total

def adam_bashforth3(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam):
	Fn=[]
	for j in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_adam[n],'y':lista_y_adam[n]}))
		n=n+1
	
	Fn_0=((23)*Fn[2])
	Fn_1=((16)*Fn[1])
	Fn_2=((5)*Fn[0])

	Fn_total=((valor_h/12)*(Fn_0-Fn_1+Fn_2))
	return Fn_total

def adam_bashforth4(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam):
	Fn=[]
	for j in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_adam[n],'y':lista_y_adam[n]}))
		n=n+1
	
	Fn_0=((55)*Fn[3])
	Fn_1=((59)*Fn[2])
	Fn_2=((37)*Fn[1])
	Fn_3=((9)*Fn[0])
	Fn_total=((valor_h/24)*(Fn_0-Fn_1+Fn_2-Fn_3))
	return Fn_total

def adam_bashforth5(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam):
	Fn=[]
	for j in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_adam[n],'y':lista_y_adam[n]}))
		n=n+1
	
	Fn_0=((1901)*Fn[4])
	Fn_1=((2774)*Fn[3])
	Fn_2=((2616)*Fn[2])
	Fn_3=((1274)*Fn[1])
	Fn_4=((251)*Fn[0])
	Fn_total=((valor_h/720)*(Fn_0-Fn_1+Fn_2-Fn_3+Fn_4))
	return Fn_total

def adam_bashforth(expressao,leitura,valor_h,valor_passos,lista_y_adam,lista_t_adam):
	n=0
	ordem=leitura-1
	valor_y=lista_y_adam[ordem]
	valor_t=lista_t_adam[ordem]

	for i in range(ordem):
		print (i , 	 lista_y_adam[i])

	for i in range(ordem,valor_passos+1):
		if ((ordem+1) == 2):
			Fn=adam_bashforth2(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 3):
			Fn=adam_bashforth3(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 4):
			Fn=adam_bashforth4(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 5):
			Fn=adam_bashforth5(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)

		n=n+1
		valor_y=valor_y + Fn
		valor_t=valor_t + valor_h
		print (i ,   valor_y)
		lista_y_adam.append(valor_y)
		lista_t_adam.append(valor_t)
	return


#MOULTON 2,3,4.....
def adam_moulton2(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
		
	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem],lista_t_moulton[ordem],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem]+valor_h
	
    
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})
	Fn=S(expressao).subs({'t':lista_t_moulton[ordem],'y':lista_y_moulton[ordem]})

	Fn_total = ((valor_h/2)*(Fn+Fn_mais_1))

	return Fn_total

def adam_moulton3(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
	Fn=[]
	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem],lista_t_moulton[ordem],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem]+valor_h
	
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	for i in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_moulton[n],'y':lista_y_moulton[n]}))
		n=n+1

	Fn_mais_1=(5*Fn_mais_1)
	Fn_0=8*Fn[1]
	Fn_1=Fn[0]

	Fn_total = ((valor_h/12)*(Fn_mais_1+Fn_0-Fn_1))

	return Fn_total

def adam_moulton4(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
	Fn=[]	

	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem],lista_t_moulton[ordem],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})


	for i in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_moulton[n],'y':lista_y_moulton[n]}))
		n=n+1

	Fn_mais_1=Fn_mais_1*9
	Fn_0=((19)*Fn[2])
	Fn_1=((5)*Fn[1])
	Fn_2=Fn[0]

	Fn_total =((valor_h/24)*(Fn_mais_1+Fn_0-Fn_1+Fn_2))
	return Fn_total

def adam_moulton(expressao,ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton):
	n=0
	ordem=ordem-1

	for k in range(ordem):
		print (k , 	 lista_y_moulton[k])

	for i in range(ordem,valor_passos+1):
		valor_y=lista_y_moulton[i]
		valor_t=lista_t_moulton[i]

		if ((ordem+1) == 2):
			Fn=adam_moulton2(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)

		if ((ordem+1) == 3):
			Fn=adam_moulton3(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		
		if ((ordem+1) == 4):
			Fn=adam_moulton4(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		
		n=n+1
		valor_y=valor_y + Fn
		valor_t=valor_t + valor_h
		print (i ,   valor_y)
		lista_y_moulton.append(valor_y)
		lista_t_moulton.append(valor_t)

	return



#FORMULA_INVERSA DE 2,3,4
def formula_inversa2(expressao,ordem,n,valor_h,lista_y_inversa,lista_t_inversa):

	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_inversa[ordem],lista_t_inversa[ordem],valor_h)	
	valor_t_mais_1=lista_t_inversa[ordem]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	Fn_mais_1=((2*valor_h)*Fn_mais_1)
	y_0=(4*(lista_y_inversa[ordem]))
	y_1=lista_y_inversa[ordem-1]

	Yn_total=((1/3)*(y_0-y_1+Fn_mais_1))

	return	Yn_total

def formula_inversa4(expressao,ordem,n,valor_h,lista_y_inversa,lista_t_inversa):

	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_inversa[ordem],lista_t_inversa[ordem],valor_h)	
	valor_t_mais_1=lista_t_inversa[ordem]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	Fn_mais_1=((12*valor_h)*Fn_mais_1)
	y_0=(48*lista_y_inversa[ordem])
	y_1=(36*lista_y_inversa[ordem-1])
	y_2=(16*lista_y_inversa[ordem-2])
	y_3=(3*lista_y_inversa[ordem-3])

	Yn_total=((1/25)*(y_0-y_1+y_2-y_3+Fn_mais_1))

	return	Yn_total

def formula_inversa(expressao,ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa):
	n=0
	ordem=ordem-1
	valor_y=lista_y_inversa[ordem]
	valor_t=lista_t_inversa[ordem]


	for k in range(ordem):
		print (k , 	 lista_y_inversa[k])


	for i in range(ordem,valor_passos+1):

		if ((ordem+1) == 2):
			Fn=formula_inversa2(expressao,ordem,n,valor_h,lista_y_inversa,lista_t_inversa)

		if ((ordem+1) == 4):
			Fn=formula_inversa2(expressao,ordem,n,valor_h,lista_y_inversa,lista_t_inversa)
	
		n=n+1
		valor_y=valor_y + Fn
		valor_t=valor_t + valor_h
		print (i ,   valor_y)
		lista_y_inversa.append(valor_y)
		lista_t_inversa.append(valor_t)

	return 










#INICIO DA ''MAIN''
#Lendo o arquivo
arquivo=open('arquivo.txt','r')
leitura=arquivo.read().split()


#Declarações Iniciais
valor_x=0
valor_y=0
valor_h=0
valor_passos=0


#Listas Usadas
lista_y_euler=[]
lista_t_euler=[]

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

lista_y_moulton=[]
lista_t_moulton=[]

lista_y_inversa=[]
lista_t_inversa=[]

#-------------------------------PASSOS SIMPLES---------------------------------#
if(leitura[0]=='euler'):
	print ('Metodo de euler')
	print ('y(',leitura[2],') =',leitura[1])
	print ('h =',leitura[3])

	lista_y_euler.append(int(leitura[1]))
	lista_t_euler.append(int(leitura[2]))
	valor_h=float(leitura[3])
	valor_passos=int(leitura[4])


	euler(leitura[5],valor_h,valor_passos,lista_y_euler,lista_t_euler,1)
	plt.plot(lista_t_euler, lista_y_euler)
	plt.show()

if(leitura[0]=='euler_inverso'):
	print ('Metodo de euler_inverso')
	print ('y(',leitura[2],') =',leitura[1])
	print ('h =',leitura[3])

	lista_y_inverso.append(int(leitura[1]))
	lista_t_inverso.append(int(leitura[2]))
	valor_h=float(leitura[3])
	valor_passos=int(leitura[4])

	euler_inverso(leitura[5],valor_h,valor_passos,lista_y_inverso,lista_t_inverso,1)
	plt.plot(lista_t_inverso, lista_y_inverso)
	plt.show()

if(leitura[0]=='euler_aprimorado'):
	
	print ('Metodo de euler_aprimorado')
	print ('y(',leitura[2],') =',leitura[1])
	print ('h =',leitura[3])

	valor_passos=int(leitura[4])
	valor_h=float(leitura[3])
	lista_y_aprimorado.append(int(leitura[1]))
	lista_t_aprimorado.append(int(leitura[2]))

	euler_aprimorado(leitura[5],valor_h,valor_passos,lista_y_aprimorado,lista_t_aprimorado,1)
	plt.plot(lista_t_aprimorado, lista_y_aprimorado)
	plt.show()

if(leitura[0]=='runge_kutta'):
	print ('Metodo de runge_kutta')
	print ('y(',leitura[2],') =',leitura[1])
	print ('h =',leitura[3])

	valor_passos=int(leitura[4])
	valor_h=float(leitura[3])
	lista_y_kutta.append(int(leitura[1]))
	lista_t_kutta.append(int(leitura[2]))

	runge_kutta(leitura[5],valor_h,valor_passos,lista_y_kutta,lista_t_kutta,1)
	#runge_kutta5(expressao,valor_y,valor_t,valor_h,valor_passos,lista_y_kutta5,lista_t_kutta5)
	plt.plot(lista_t_kutta, lista_y_kutta)
	plt.show()


#-------------------------------ADAM_BASHFORTH-----------------------------#
if(leitura[0]=='adam_bashforth'):	
		t3=float(leitura[len(leitura)-5])
		valor_h=float(leitura[len(leitura)-4])
		
		for i in range(1,len(leitura)-5):
			y3=float(leitura[i])
			lista_y_adam.append(y3)
			lista_t_adam.append(t3)
			t3=t3+valor_h

		valor_passos=int(leitura[len(leitura)-3])
		ordem=int(leitura[len(leitura)-1])
		adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)
		plt.plot(lista_t_adam, lista_y_adam)
		plt.show()

if(leitura[0]=='adam_bashforth_by_euler'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_adam.append(int(leitura[1]))
	lista_t_adam.append(int(leitura[2]))

	euler(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,0)
	adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)

	plt.plot(lista_t_adam, lista_y_adam)
	plt.show()

if(leitura[0]=='adam_bashforth_by_euler_inverso'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_adam.append(int(leitura[1]))
	lista_t_adam.append(int(leitura[2]))

	euler_inverso(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,0)
	adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)

	plt.plot(lista_t_adam, lista_y_adam)
	plt.show()

if(leitura[0]=='adam_bashforth_by_euler_aprimorado'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_adam.append(int(leitura[1]))
	lista_t_adam.append(int(leitura[2]))

	euler_aprimorado(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,0)
	adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)

	plt.plot(lista_t_adam, lista_y_adam)
	plt.show()

if(leitura[0]=='adam_bashforth_by_runge_kutta'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_adam.append(int(leitura[1]))
	lista_t_adam.append(int(leitura[2]))

	runge_kutta(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,0)
	adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)

	plt.plot(lista_t_adam, lista_y_adam)
	plt.show()

#------------------------------ADAM MOULTON---------------------------------#
if(leitura[0]=='adam_multon'):	
	t3=float(leitura[len(leitura)-5])
	valor_h=float(leitura[len(leitura)-4])
		
	for i in range(1,len(leitura)-5):
		y3=float(leitura[i])
		lista_y_moulton.append(y3)
		lista_t_moulton.append(t3)
		t3=t3+valor_h

	valor_passos=int(leitura[len(leitura)-3])
	ordem=int(leitura[len(leitura)-1])
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)
	plt.plot(lista_t_moulton, lista_y_moulton)
	plt.show()

if(leitura[0]=='adam_multon_by_euler'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_moulton.append(int(leitura[1]))
	lista_t_moulton.append(int(leitura[2]))

	euler(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)

	plt.plot(lista_t_moulton, lista_y_moulton)
	plt.show()

if(leitura[0]=='adam_multon_by_euler_inverso'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_moulton.append(int(leitura[1]))
	lista_t_moulton.append(int(leitura[2]))

	euler_inverso(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)

	plt.plot(lista_t_moulton, lista_y_moulton)
	plt.show()

if(leitura[0]=='adam_multon_by_euler_aprimorado'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_moulton.append(int(leitura[1]))
	lista_t_moulton.append(int(leitura[2]))

	euler_aprimorado(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)

	plt.plot(lista_t_moulton, lista_y_moulton)
	plt.show()

if(leitura[0]=='adam_multon_by_runge_kutta'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_moulton.append(int(leitura[1]))
	lista_t_moulton.append(int(leitura[2]))

	runge_kutta(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)

	plt.plot(lista_t_moulton, lista_y_moulton)
	plt.show()

#----------------------------FORMULA INVERSA----------------------------------#
if(leitura[0]=='formula_inversa'):
		t3=float(leitura[len(leitura)-5])
		valor_h=float(leitura[len(leitura)-4])
		
		for i in range(1,len(leitura)-5):
			y3=float(leitura[i])
			lista_y_inversa.append(y3)
			lista_t_inversa.append(t3)
			t3=t3+valor_h

		valor_passos=int(leitura[len(leitura)-3])
		ordem=int(leitura[len(leitura)-1])
		formula_inversa(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)
		plt.plot(lista_t_inversa, lista_y_inversa)
		plt.show()

if(leitura[0]=='formula_inversa_by_euler'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_inversa.append(int(leitura[1]))
	lista_t_inversa.append(int(leitura[2]))

	euler(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)

	plt.plot(lista_t_inversa, lista_y_inversa)
	plt.show()

if(leitura[0]=='formula_inversa_by_euler_inverso'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_inversa.append(int(leitura[1]))
	lista_t_inversa.append(int(leitura[2]))

	euler_inverso(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)

	plt.plot(lista_t_inversa, lista_y_inversa)
	plt.show()

if(leitura[0]=='formula_inversa_by_euler_aprimorado'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_inversa.append(int(leitura[1]))
	lista_t_inversa.append(int(leitura[2]))

	euler_aprimorado(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)

	plt.plot(lista_t_inversa, lista_y_inversa)
	plt.show()

if(leitura[0]=='formula_inversa_by_runge_kutta'):	
	ordem=int(leitura[len(leitura)-1])
	valor_passos=int(leitura[len(leitura)-3])
	valor_h=float(leitura[3])
	lista_y_inversa.append(int(leitura[1]))
	lista_t_inversa.append(int(leitura[2]))

	runge_kutta(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,0)
	adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)

	plt.plot(lista_t_inversa, lista_y_inversa)
	plt.show()



arquivo.close()