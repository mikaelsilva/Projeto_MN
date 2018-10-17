from sympy import symbols, Function, solve,Eq,S
from matplotlib import pyplot as plt
import sympy 



#---------------------------------------------------------------------------------#
#OK
def euler(expressao,valor_h,valor_passos,lista_y_euler,lista_t_euler,euler):
	
	valor_y=lista_y_euler[0]
	valor_t=lista_t_euler[0]

	for i in range(valor_passos-euler):
		Fn=S(expressao).subs({'t':lista_t_euler[i],'y':lista_y_euler[i]})
		valor_y=valor_y+valor_h*Fn
		valor_t=valor_t+valor_h

		lista_y_euler.append(valor_y)
		lista_t_euler.append(valor_t)	
	
	return 


def euler_n_mais_1(expressao,valor_y2,valor_t2,valor_h2):

	Fn2=S(expressao).subs({'t':valor_t2,'y':valor_y2})
	valor_y2=valor_y2+valor_h2*Fn2
	
	return valor_y2

#OK
def euler_inverso(expressao,valor_h,valor_passos,lista_y_inverso,lista_t_inverso,inverso):
	
	valor_y=lista_y_inverso[0]
	valor_t=lista_t_inverso[0]

	
	for i in range(valor_passos-inverso):
		
		valor_y_mais_1=euler_n_mais_1(expressao,valor_y,valor_t,valor_h) #Yn+1 para Fn+1
		valor_t_mais_1=valor_t+valor_h #Tn+1

		Fn=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1}) #Fn+1
		valor_y=valor_y+valor_h*Fn
		valor_t=valor_t+valor_h

		lista_y_inverso.append(valor_y)
		lista_t_inverso.append(valor_t)

	return 

#OK 
def euler_aprimorado(expressao,valor_h,valor_passos,lista_y_aprimorado,lista_t_aprimorado,aprimorado):

	valor_y=lista_y_aprimorado[0]
	valor_t=lista_t_aprimorado[0]
	
	for i in range(valor_passos-aprimorado):

		#Fn_mais_1
		Fn_mais_1=S(expressao).subs({'t':valor_t,'y':valor_y})
		valor_y_mais_1=valor_y+valor_h*Fn_mais_1
		valor_t_mais_1=valor_t+valor_h
		Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

		#Fn
		Fn=S(expressao).subs({'t':lista_t_aprimorado[i],'y':lista_y_aprimorado[i]})
		Fn_total=Fn+Fn_mais_1

		#Yn+1
		valor_y=valor_y+(valor_h/2)*Fn_total
		valor_t=valor_t+valor_h

		#Guardando na Lista
		lista_y_aprimorado.append(valor_y)
		lista_t_aprimorado.append(valor_t)

	return

#OK
def runge_kutta(expressao,valor_h,valor_passos,lista_y_kutta,lista_t_kutta,kutta):


	valor_y=lista_y_kutta[0]
	valor_t=lista_t_kutta[0]

	for i in range(valor_passos-kutta):

		#Coeficientes para encontrar o Fn
		K_1=S(expressao).subs({'t':valor_t,'y':valor_y})
		K_2=S(expressao).subs({'t': valor_t+(0.5)*valor_h ,'y':valor_y+(0.5*valor_h*K_1)})
		K_3=S(expressao).subs({'t': valor_t+(0.5)*valor_h ,'y':valor_y+(0.5*valor_h*K_2)})
		K_4=S(expressao).subs({'t': valor_t+valor_h ,'y':valor_y+valor_h*K_3})

		#Calcula o Fn e o Yn+1 e Tn+1
		Fn=(K_1+ 2*K_2 + 2*K_3 + K_4)
		valor_y=valor_y+(valor_h/6)*Fn
		valor_t=valor_t+valor_h

		#Guarda os valores
		lista_y_kutta.append(valor_y)
		lista_t_kutta.append(valor_t)

	return 

#VERIFICAR
def runge_kutta5(expressao,valor_h,valor_passos,valor_y_kutta5,valor_t_kutta5,kutta5):
		
	valor_t=lista_t_kutta5[0]
	valor_y=lista_y_kutta5[0]
	
	for next in range(valor_passos-kutta5):

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

		lista_y_kutta5.append(valor_y)
		lista_t_kutta5.append(valor_t)
	return 

#---------------------------------------------------------------------------------#
#BASHFORTH 2,3,4,5,6,7,8
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

def adam_bashforth6(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam):
	Fn=[]
	for j in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_adam[n],'y':lista_y_adam[n]}))
		n=n+1
	
	Fn_0=((4277)*Fn[5])
	Fn_1=((7923)*Fn[4])
	Fn_2=((9982)*Fn[3])
	Fn_3=((7298)*Fn[2])
	Fn_4=((2877)*Fn[1])
	Fn_5=((475)*(Fn[0]))
	Fn_total=((valor_h/1440)*(Fn_0-Fn_1+Fn_2-Fn_3+Fn_4-Fn_5))

	return Fn_total

def adam_bashforth7(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam):
	Fn=[]
	for j in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_adam[n],'y':lista_y_adam[n]}))
		n=n+1
	
	Fn_0=((198721)*Fn[6])
	Fn_1=((447288)*Fn[5])
	Fn_2=((705549)*Fn[4])
	Fn_3=((688256)*Fn[3])
	Fn_4=((407139)*Fn[2])
	Fn_5=((134472)*(Fn[1]))
	Fn_6=((19087)*(Fn[0]))
	Fn_total=((valor_h/60480)*(Fn_0-Fn_1+Fn_2-Fn_3+Fn_4-Fn_5+Fn_6))
	return Fn_total

def adam_bashforth8(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam):
	Fn=[]
	for j in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_adam[n],'y':lista_y_adam[n]}))
		n=n+1
	
	Fn_0=((434241)*Fn[7])
	Fn_1=((1152169)*Fn[6])
	Fn_2=((2183877)*Fn[5])
	Fn_3=((2664477)*Fn[4])
	Fn_4=((2102243)*Fn[3])
	Fn_5=((1041723)*(Fn[2]))
	Fn_6=((295767)*(Fn[1]))
	Fn_7=((36799)*(Fn[0]))
	Fn_total=((valor_h/120960)*(Fn_0-Fn_1+Fn_2-Fn_3+Fn_4-Fn_5+Fn_6-Fn_7))
	return Fn_total

def adam_bashforth(expressao,ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam):
	#n=0 durante ordem do adam, n+1 durante ordem do adam.........
	n=0	

	ordem=ordem-1
	valor_y=lista_y_adam[ordem]
	valor_t=lista_t_adam[ordem]

	for i in range(ordem+1,valor_passos+1):
		if ((ordem+1) == 2):
			Fn=adam_bashforth2(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 3):
			Fn=adam_bashforth3(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 4):
			Fn=adam_bashforth4(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 5):
			Fn=adam_bashforth5(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 6):
			Fn=adam_bashforth6(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 7):
			Fn=adam_bashforth7(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)
		if ((ordem+1) == 8):
			Fn=adam_bashforth8(expressao,ordem,n,valor_h,lista_y_adam,lista_t_adam)

		n=n+1
		valor_y=valor_y + Fn
		valor_t=valor_t + valor_h

		lista_y_adam.append(valor_y)
		lista_t_adam.append(valor_t)
	return


#---------------------------------------------------------------------------------#
#MOULTON 2,3,4.....
def adam_moulton2(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
		
	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem],lista_t_moulton[ordem],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem]+valor_h
	
    #Fn+1
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})
	
	#Fn
	Fn=S(expressao).subs({'t':lista_t_moulton[ordem],'y':lista_y_moulton[ordem]})

	#Fn final
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

def adam_moulton5(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
	Fn=[]	

	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem],lista_t_moulton[ordem],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})


	for i in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_moulton[n],'y':lista_y_moulton[n]}))
		n=n+1

	Fn_mais_1=Fn_mais_1*251
	Fn_0=((646)*Fn[3])
	Fn_1=((264)*Fn[2])
	Fn_2=((106)*Fn[1])
	Fn_3=((19)*Fn[0])
	Fn_total =((valor_h/720)*(Fn_mais_1+Fn_0-Fn_1+Fn_2-Fn_3))
	return Fn_total

def adam_moulton6(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
	Fn=[]	

	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem+n],lista_t_moulton[ordem+n],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem+n]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})


	for i in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_moulton[n],'y':lista_y_moulton[n]}))
		n=n+1

	Fn_mais_1=Fn_mais_1*475
	Fn_0=((1427)*Fn[4])
	Fn_1=((798)*Fn[3])
	Fn_2=((482)*Fn[2])
	Fn_3=((173)*Fn[1])
	Fn_4=((27)*Fn[0])

	Fn_total =((valor_h/1440)*(Fn_mais_1+Fn_0-Fn_1+Fn_2-Fn_3+Fn_4))


	return Fn_total

def adam_moulton7(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
	Fn=[]	

	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem],lista_t_moulton[ordem],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})


	for i in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_moulton[n],'y':lista_y_moulton[n]}))
		n=n+1

	Fn_mais_1=Fn_mais_1*19087
	Fn_0=((65112)*Fn[5])
	Fn_1=((46461)*Fn[4])
	Fn_2=((37504)*Fn[3])
	Fn_3=((20211)*Fn[2])
	Fn_4=((6312)*Fn[1])
	Fn_5=((863)*Fn[0])
	Fn_total =((valor_h/60480)*(Fn_mais_1+Fn_0-Fn_1+Fn_2-Fn_3+Fn_4-Fn_5))
	return Fn_total

def adam_moulton8(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton):
	Fn=[]	

	#Calcula o Y+1 por euler	
	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_moulton[ordem],lista_t_moulton[ordem],valor_h)	
	valor_t_mais_1=lista_t_moulton[ordem]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})


	for i in range(ordem+1):
		Fn.append(S(expressao).subs({'t':lista_t_moulton[n],'y':lista_y_moulton[n]}))
		n=n+1

	Fn_mais_1=Fn_mais_1*36799
	Fn_0=((139849)*Fn[6])
	Fn_1=((121797)*Fn[5])
	Fn_2=((123133)*Fn[4])
	Fn_3=((88547)*Fn[3])
	Fn_4=((41499)*Fn[2])
	Fn_5=((11351)*Fn[1])
	Fn_6=((1375)*Fn[0])
	Fn_total =((valor_h/120960)*(Fn_mais_1+Fn_0-Fn_1+Fn_2-Fn_3+Fn_4-Fn_6+Fn_6))
	return Fn_total

def adam_moulton(expressao,ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton):
	n=0
	ordem=ordem-1
	
	
	for i in range(ordem+1,valor_passos+1):
		valor_y=lista_y_moulton[i-1]
		valor_t=lista_t_moulton[i-1]

		if ((ordem+1) == 2):
			Fn=adam_moulton2(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		if ((ordem+1) == 3):
			Fn=adam_moulton3(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		if ((ordem+1) == 4):
			Fn=adam_moulton4(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		if ((ordem+1) == 5):
			Fn=adam_moulton5(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		if ((ordem+1) == 6):
			Fn=adam_moulton6(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		if ((ordem+1) == 7):
			Fn=adam_moulton7(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)
		if ((ordem+1) == 8):
			Fn=adam_moulton8(expressao,ordem,n,valor_h,lista_y_moulton,lista_t_moulton)

		n=n+1
		valor_y=valor_y + Fn
		valor_t=valor_t + valor_h

		lista_y_moulton.append(valor_y)
		lista_t_moulton.append(valor_t)

	return



#---------------------------------------------------------------------------------#
#FORMULA_INVERSA DE 2,3,4
def formula_inversa2(expressao,n,valor_h,lista_y_inversa,lista_t_inversa):

	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_inversa[n],lista_t_inversa[n],valor_h)	
	valor_t_mais_1=lista_t_inversa[n]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	Fn_mais_1=((2*valor_h)*Fn_mais_1)
	y_0=(4*(lista_y_inversa[n]))
	y_1=lista_y_inversa[n-1]

	Yn_total=((1/3)*(y_0-y_1+Fn_mais_1))

	return	Yn_total

def formula_inversa3(expressao,n,valor_h,lista_y_inversa,lista_t_inversa):

	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_inversa[n],lista_t_inversa[n],valor_h)	
	valor_t_mais_1=lista_t_inversa[n]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	Fn_mais_1=((6*valor_h)*Fn_mais_1)
	y_0=(18*lista_y_inversa[n])
	y_1=(9*lista_y_inversa[n-1])
	y_2=(2*lista_y_inversa[n-2])

	Yn_total=((1/11)*(y_0-y_1+y_2+Fn_mais_1))

	return	Yn_total

def formula_inversa4(expressao,n,valor_h,lista_y_inversa,lista_t_inversa):

	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_inversa[n],lista_t_inversa[n],valor_h)	
	valor_t_mais_1=lista_t_inversa[n]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	Fn_mais_1=((12*valor_h)*Fn_mais_1)
	y_0=(48*lista_y_inversa[n])
	y_1=(36*lista_y_inversa[n-1])
	y_2=(16*lista_y_inversa[n-2])
	y_3=(3*lista_y_inversa[n-3])

	Yn_total=((1/25)*(y_0-y_1+y_2-y_3+Fn_mais_1))

	return	Yn_total

def formula_inversa5(expressao,n,valor_h,lista_y_inversa,lista_t_inversa):

	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_inversa[n],lista_t_inversa[n],valor_h)	
	valor_t_mais_1=lista_t_inversa[n]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	Fn_mais_1=((60*valor_h)*Fn_mais_1)
	y_0=(300*lista_y_inversa[n])
	y_1=(300*lista_y_inversa[n-1])
	y_2=(200*lista_y_inversa[n-2])
	y_3=(75*lista_y_inversa[n-3])
	y_4=(12*lista_y_inversa[n-4])

	Yn_total=((1/137)*(y_0-y_1+y_2-y_3+y_4+Fn_mais_1))

	return	Yn_total

def formula_inversa6(expressao,n,valor_h,lista_y_inversa,lista_t_inversa):

	valor_y_mais_1=euler_n_mais_1(expressao,lista_y_inversa[n],lista_t_inversa[n],valor_h)	
	valor_t_mais_1=lista_t_inversa[n]+valor_h
	Fn_mais_1=S(expressao).subs({'t':valor_t_mais_1,'y':valor_y_mais_1})

	Fn_mais_1=((60*valor_h*Fn_mais_1))
	y_0=(360*lista_y_inversa[n])
	y_1=(450*lista_y_inversa[n-1])
	y_2=(400*lista_y_inversa[n-2])
	y_3=(225*lista_y_inversa[n-3])
	y_4=(72*lista_y_inversa[n-4])
	y_5=(10*lista_y_inversa[n-5])

	Yn_total=((y_0-y_1+y_2-y_3+y_4-y_5+Fn_mais_1)/147)

	return	Yn_total

def formula_inversa(expressao,ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa):
	
	ordem=ordem-1
	n=ordem
	valor_y=lista_y_inversa[ordem]
	valor_t=lista_t_inversa[ordem]
	

	for i in range(ordem+1,valor_passos+1):
		
		if ((ordem+1) == 2):
			Fn=formula_inversa2(expressao,n,valor_h,lista_y_inversa,lista_t_inversa)
		if ((ordem+1) == 3):
			Fn=formula_inversa3(expressao,n,valor_h,lista_y_inversa,lista_t_inversa)
		if ((ordem+1) == 4):
			Fn=formula_inversa4(expressao,n,valor_h,lista_y_inversa,lista_t_inversa)
		if ((ordem+1) == 5):
			Fn=formula_inversa5(expressao,n,valor_h,lista_y_inversa,lista_t_inversa)
		if ((ordem+1) == 6):
			Fn=formula_inversa6(expressao,n,valor_h,lista_y_inversa,lista_t_inversa)
	

		n=n+1
		valor_y=Fn 	   #Errado com fator de 0,0999...
		valor_t=valor_t + valor_h

		lista_y_inversa.append(valor_y)
		lista_t_inversa.append(valor_t)


	return 





#INICIO DA ''MAIN''


#Declaracoes Iniciais
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







salvar=open('salvar.txt','w')
arquivo=open('arquivo.txt','r')
with open('arquivo.txt') as final:
	for i in final:
		leitura=arquivo.readline().split()
		

#-------------------------------PASSOS SIMPLES---------------------------------#
		if(leitura[0]=='euler'):
			salvar.write('Metodo de euler\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[3])+ '\n')



			lista_y_euler.append(float(leitura[1]))
			lista_t_euler.append(float(leitura[2]))
			valor_h=float(leitura[3])
			valor_passos=int(leitura[4])


			euler(salvar,leitura[5],valor_h,valor_passos,lista_y_euler,lista_t_euler,0)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_euler[i]) + '\n')


			plt.plot(lista_t_euler, lista_y_euler)
			plt.show()
		

			lista_t_euler.clear()
			lista_y_euler.clear()

		if(leitura[0]=='euler_inverso'):
			salvar.write('Metodo de Euler_Inverso\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[3])+ '\n')

			lista_y_inverso.append(float(leitura[1]))
			lista_t_inverso.append(float(leitura[2]))
			valor_h=float(leitura[3])
			valor_passos=int(leitura[4])

			euler_inverso(leitura[5],valor_h,valor_passos,lista_y_inverso,lista_t_inverso,0)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_inverso[i]) + '\n')

			plt.plot(lista_t_inverso, lista_y_inverso)
			plt.show()
			
			lista_t_inverso.clear()
			lista_y_inverso.clear()

		if(leitura[0]=='euler_aprimorado'):
			salvar.write('Metodo de Euler_Aprimorado\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[3])+ '\n')

			valor_passos=int(leitura[4])
			valor_h=float(leitura[3])
			lista_y_aprimorado.append(float(leitura[1]))
			lista_t_aprimorado.append(float(leitura[2]))

			euler_aprimorado(leitura[5],valor_h,valor_passos,lista_y_aprimorado,lista_t_aprimorado,0)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_aprimorado[i]) + '\n')

			plt.plot(lista_t_aprimorado, lista_y_aprimorado)
			plt.show()


			lista_t_aprimorado.clear()
			lista_y_aprimorado.clear()

		if(leitura[0]=='runge_kutta'):
			salvar.write('Metodo de Runge_Kutta\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[3])+ '\n')

			valor_passos=int(leitura[4])
			valor_h=float(leitura[3])
			lista_y_kutta.append(float(leitura[1]))
			lista_t_kutta.append(float(leitura[2]))

			runge_kutta(leitura[5],valor_h,valor_passos,lista_y_kutta,lista_t_kutta,0)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_kutta[i]) + '\n')

			plt.plot(lista_t_kutta, lista_y_kutta)
			plt.show()

			lista_t_kutta.clear()
			lista_y_kutta.clear()

		if(leitura[0]=='runge_kutta5'):
			salvar.write('Metodo de Runge_Kutta5\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[3])+ '\n')

			valor_passos=int(leitura[4])
			valor_h=float(leitura[3])
			lista_y_kutta5.append(float(leitura[1]))
			lista_t_kutta5.append(float(leitura[2]))

			runge_kutta5(leitura[5],valor_h,valor_passos,lista_y_kutta5,lista_t_kutta5,0)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_kutta5[i]) + '\n')


			plt.plot(lista_t_kutta5, lista_y_kutta5)
			plt.show()

			lista_t_kutta5.clear()
			lista_y_kutta5.clear()
		
#-------------------------------ADAM_BASHFORTH-----------------------------#
		if(leitura[0]=='adam_bashforth'):	
			salvar.write('Metodo de Adam_Bashforth\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

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
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_adam[i]) + '\n')

			plt.plot(lista_t_adam, lista_y_adam)
			plt.show()
		
		if(leitura[0]=='adam_bashforth_by_euler'):	
			salvar.write('Metodo de Adam_Bashforth_By_Euler\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_adam.append(float(leitura[1]))
			lista_t_adam.append(float(leitura[2]))

			euler(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,1)
			adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_adam[i]) + '\n')


			plt.plot(lista_t_adam, lista_y_adam)
			plt.show()

		if(leitura[0]=='adam_bashforth_by_euler_inverso'):	
			salvar.write('Metodo de Adam_Bashforth_By_Euler_Inverso\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_adam.append(float(leitura[1]))
			lista_t_adam.append(float(leitura[2]))

			euler_inverso(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,1)
			adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_adam[i]) + '\n')

			plt.plot(lista_t_adam, lista_y_adam)
			plt.show()

		if(leitura[0]=='adam_bashforth_by_euler_aprimorado'):	
			salvar.write('Metodo de Adam_Bashforth_By_Euler_Aprimorado\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_adam.append(float(leitura[1]))
			lista_t_adam.append(float(leitura[2]))

			euler_aprimorado(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,1)
			adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_adam[i]) + '\n')

			plt.plot(lista_t_adam, lista_y_adam)
			plt.show()

		if(leitura[0]=='adam_bashforth_by_runge_kutta'):	
			salvar.write('Metodo de Adam_Bashforth_By_Runge_Kutta\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_adam.append(float(leitura[1]))
			lista_t_adam.append(float(leitura[2]))

			runge_kutta(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,1)
			adam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_adam[i]) + '\n')
			plt.plot(lista_t_adam, lista_y_adam)
			plt.show()


		#CORRIGIRRRR NÃO ESQUECERRR
		if(leitura[0]=='adam_bashforth_by_runge_kutta5'):	
			salvar.write('Metodo de Adam_Bashforth_By_Runge_Kutta5\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_adam.append(float(leitura[1]))
			lista_t_adam.append(float(leitura[2]))

			#runge_kutta5(leitura[len(leitura)-2],valor_h,ordem,lista_y_adam,lista_t_adam,1)
			#dam_bashforth(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_adam,lista_t_adam)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_adam[i]) + '\n')
			plt.plot(lista_t_adam, lista_y_adam)
			plt.show()

		lista_t_adam.clear()
		lista_y_adam.clear()
#------------------------------ADAM MOULTON---------------------------------#
		if(leitura[0]=='adam_multon'):	
			salvar.write('Metodo de Adam_Multon\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

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
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_moulton[i]) + '\n')

			plt.plot(lista_t_moulton, lista_y_moulton)
			plt.show()

		if(leitura[0]=='adam_multon_by_euler'):	
			salvar.write('Metodo de Adam_Multon_By_Euler\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_moulton.append(float(leitura[1]))
			lista_t_moulton.append(float(leitura[2]))

			euler(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,1)
			adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_moulton[i]) + '\n')

			plt.plot(lista_t_moulton, lista_y_moulton)
			plt.show()

		if(leitura[0]=='adam_multon_by_euler_inverso'):	
			print ('Metodo de adam multon by euler inverso')
			print ('y(',leitura[2],') =',leitura[1])
			print ('h =',leitura[len(leitura)-4])
			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_moulton.append(float(leitura[1]))
			lista_t_moulton.append(float(leitura[2]))

			euler_inverso(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,0)
			adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)

			plt.plot(lista_t_moulton, lista_y_moulton)
			plt.show()

		if(leitura[0]=='adam_multon_by_euler_aprimorado'):	
			print ('Metodo de adam multon by euler aprimorado')
			print ('y(',leitura[2],') =',leitura[1])
			print ('h =',leitura[len(leitura)-4])
			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_moulton.append(float(leitura[1]))
			lista_t_moulton.append(float(leitura[2]))

			euler_aprimorado(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,0)
			adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)

			plt.plot(lista_t_moulton, lista_y_moulton)
			plt.show()

		if(leitura[0]=='adam_multon_by_runge_kutta'):	
			print ('Metodo de adam multon by runge kutta')
			print ('y(',leitura[2],') =',leitura[1])
			print ('h =',leitura[len(leitura)-4])
			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_moulton.append(float(leitura[1]))
			lista_t_moulton.append(float(leitura[2]))

			runge_kutta(leitura[len(leitura)-2],valor_h,ordem,lista_y_moulton,lista_t_moulton,0)
			adam_moulton(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_moulton,lista_t_moulton)

			plt.plot(lista_t_moulton, lista_y_moulton)
			plt.show()

		lista_t_moulton.clear()
		lista_y_moulton.clear()
#----------------------------FORMULA INVERSA----------------------------------#
		if(leitura[0]=='formula_inversa'):
			salvar.write('Metodo de Formula_Inversa\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

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
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_inversa[i]) + '\n')

			plt.plot(lista_t_inversa, lista_y_inversa)
			plt.show()

		if(leitura[0]=='formula_inversa_by_euler'):	
			salvar.write('Metodo de Formula_Inversa_By_Euler\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_inversa.append(float(leitura[1]))
			lista_t_inversa.append(float(leitura[2]))

			euler(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,1)
			formula_inversa(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_inversa[i]) + '\n')

			plt.plot(lista_t_inversa, lista_y_inversa)
			plt.show()

		if(leitura[0]=='formula_inversa_by_euler_inverso'):	
			salvar.write('Metodo de Formula_Inversa_By_Euler_Inverso\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n')

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_inversa.append(float(leitura[1]))
			lista_t_inversa.append(float(leitura[2]))

			euler_inverso(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,1)
			formula_inversa(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_inversa[i]) + '\n')
			
			plt.plot(lista_t_inversa, lista_y_inversa)
			plt.show()

		if(leitura[0]=='formula_inversa_by_euler_aprimorado'):	
			salvar.write('Metodo de Formula_Inversa_By_Euler_Aprimorado\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n') 

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_inversa.append(float(leitura[1]))
			lista_t_inversa.append(float(leitura[2]))

			euler_aprimorado(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,1)
			formula_inversa(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_inversa[i]) + '\n')

			plt.plot(lista_t_inversa, lista_y_inversa)
			plt.show()

		if(leitura[0]=='formula_inversa_by_runge_kutta'):	
			salvar.write('Metodo de Formula_Inversa_By_Runge_Kutta\n')
			salvar.write('y(' + str(leitura[2]) +') = ' + str(leitura[1]) + '\n')
			salvar.write('h =  '+ str(leitura[len(leitura)-4])+ '\n') 

			ordem=int(leitura[len(leitura)-1])
			valor_passos=int(leitura[len(leitura)-3])
			valor_h=float(leitura[3])
			lista_y_inversa.append(float(leitura[1]))
			lista_t_inversa.append(float(leitura[2]))

			runge_kutta(leitura[len(leitura)-2],valor_h,ordem,lista_y_inversa,lista_t_inversa,1)
			formula_inversa(leitura[len(leitura)-2],ordem,valor_h,valor_passos,lista_y_inversa,lista_t_inversa)
			for i in range(21):
				salvar.write(str(i) + ' ' + str(lista_y_inversa[i]) + '\n')

			plt.plot(lista_t_inversa, lista_y_inversa)
			plt.show()

		lista_t_inversa.clear()
		lista_y_inversa.clear()

arquivo.close()
salvar.close()

