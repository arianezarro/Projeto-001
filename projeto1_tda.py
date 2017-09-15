#importando os pacotes necessários.

import numpy as np #manipulação de matrizes.
import matplotlib.pyplot as plt #plot das imagens. 
from astropy.io import fits #leitura dos arquivos fits.
import glob #listar arquivos de diretórios.

#1)BIAS

def bias():#função de cálculo do bias
    arraymbias =np.array([])  #array vazio para ser usado fora do laço for.
    master_bias = [] #variável vazia que ainda será utilizada.
    for i in glob.glob('/home/ariane/fitsfile/bias.B.*.fits'): #glob-leitura dos arquivos padronizados. Nesse caso, só serão lidos os arquivos com o nome 'bias.B.*.fits' sendo * o número da imagem. fitsfile é o diretório em que estão os arquivos fits.
        b_img, hdr = fits.getdata(i, header=True) #lendo as imagens bias e seus headers. 
        b_img = np.array(b_img, dtype='Float64') #só para inserir o 'float64' mesmo.
        master_bias.append(b_img) #juntar as imagens bias do filtro B.
    master_bias = np.array(master_bias) #transformar em array.
    arraymbias2 = np.vstack([arraymbias, master_bias]) if arraymbias.size else master_bias # colocar o array 'master_bias' dentro do array vazio 'arraymbias'que está fora do laço for. Fizemos isso com a função 'vstack' do numpy.
    arraymbias2 = np.array(arraymbias2, dtype='Float64') #só para inserir o 'float64' mesmo.
    b_mediana = np.median(arraymbias2, axis=0) #calcular a mediana d imagem de bias.
    b_mediana = np.array(b_mediana, dtype='Float64') #só para inserir o 'float64' mesmo.
    plt.figure() #plotar o bias
    plt.imshow(b_mediana,vmin=np.mean(b_mediana)-2.5*np.std(b_mediana),vmax=np.mean(b_mediana)+2.5*np.std(b_mediana), cmap=plt.cm.gray,origin='lower')     
    plt.colorbar()
    plt.show()  
    return(b_mediana) #retornar a matriz final do bias.

#2)FLATFIELD

def flat_bias(bias): #função de cálculo do flat final com o bias retirado. Por isso o input será bias.
    bias1 = np.array(bias(), dtype = 'Float64') 
    arraymflat =np.array([])  #array vazio para ser usado fora do laço for.
    master_flat = [] #variável vazia que ainda será utilizada.
    for i in glob.glob('/home/ariane/fitsfile/flat.B.*.fits'): #glob-leitura dos arquivos padronizados. Nesse caso, só serão lidos os arquivos com o nome 'flat.B.*.fits' sendo * o número da imagem. fitsfile é o diretório em que estão os arquivos fits.
        f_img, hdr = fits.getdata(i, header=True) #lendo as imagens flat e seus headers. 
        f_img = np.array(f_img, dtype='Float64') #só para inserir o 'float64' mesmo.
# primeiro devemos subtrair o bias do flat para depois fazer o master flat e a normalização.
        flat_bias = f_img - bias1 #subtraindo o bias final de cada imagem de flat.
        flat_bias = np.array(flat_bias, dtype='Float64')  
# agora podemos normalizar,ou seja, dividir o flat(com bias subtraído) pela média do flat(com bias subtraído) e depois juntar as imagens flat num master flat.
        flat_bias_media = np.mean(flat_bias)
        flat_norm = flat_bias/flat_bias_media
        master_flat.append(flat_norm)
    master_flat = np.array(master_flat) #transformar em array.
    arraymflat2 = np.vstack([arraymflat, master_flat]) if arraymflat.size else master_flat # colocar o array 'master_flat' dentro do array vazio 'arraymflat'que está fora do laço for. Fizemos isso com a função 'vstack' do numpy.
    arraymflat2 = np.array(arraymflat2, dtype='Float64')
    f_mediana = np.median(arraymflat2, axis=0)
    plt.figure() #plotar o bias
    plt.imshow(f_mediana,vmin=np.mean(f_mediana)-2.5*np.std(f_mediana),vmax=np.mean(f_mediana)+2.5*np.std(f_mediana), cmap=plt.cm.gray,origin='lower')     
    plt.colorbar()
    plt.show()  
    return(f_mediana) #retornar a matriz final do bias.

#3)BACKGROUND

def sky(bias, flat_bias):#cáculo da imagem do céu a partir da imagem corrigida pelo bias e o flat.
    bias1 = np.array(bias(), dtype = 'Float64') 
    flat_bias1 = np.array(flat_bias(bias), dtype = 'Float64')
    sky_img, hdr = fits.getdata(input('Imagem de ciência para retirar o sky image: '), header = input('Header da imagem: True ou False? ')) #usuário escolherá a imagem e indicará se o header está presente ou não.
    sky_img = np.array(sky_img ,dtype = 'Float64')
    sky_bias = sky_img - bias1 #tirar o bias da imagem de ciência.
    sky_bias_flat = sky_bias/flat_bias1
        
#agora temos que perguntar ao usuário que valores de x e y ele deseja 'recortar' da imagem de ciência para fazer a imagem do céu.
    size = sky_img.shape #dimensão do array.
    print('Escolha os valores de x e y entre 0 e %s: '%str(size[1]))
    xi = int(input('xi= '))
    xf = int(input('xf= '))
    yi = int(input('yi= '))
    yf = int(input('yf= '))
    sky_bias_flat = sky_bias_flat[yi:yf, xi:xf]
    
#depois de dados os valores, calcular a função de poisson correspondente.
    sky_image = np.random.poisson(np.mean(sky_bias_flat), size) #size: int or tuple of ints.
    plt.figure() #plotar a imagem do céu.
    plt.imshow(sky_image,vmin=np.mean(sky_image)-2.5*np.std(sky_image),vmax=np.mean(sky_image)+2.5*np.std(sky_image), cmap=plt.cm.gray,origin='lower')     
    plt.colorbar()
    plt.show()  
    return(sky_image) #retornar a matriz final da imagem do céu.


#4)CORREÇÃO DAS IMAGENS DE CIÊNCIA

def img_final(bias, flat_bias, sky): #cálculo das imagens de ciência corrigidas.
    bias1 = np.array(bias(), dtype = 'Float64') 
    flat_bias1 = np.array(flat_bias(bias), dtype = 'Float64')
    sky_bias_flat1 = np.array(sky(bias,flat_bias), dtype = 'Float64')
    for i in glob.glob('/home/ariane/fitsfile/xo2b.*.fits'):
        img_fin, hdr = fits.getdata(i, header=True)
        img_fin = np.array(img_fin, dtype='Float64')
        hdr = fits.getheader(i) #salvar imagem.  
        outfile = '_bfs_corrigida.fits'

        b_img_fin = img_fin - bias1 #corrigir o bias da imagem de ciência.
        f_b_img_fin = b_img_fin/flat_bias1 #corrigir o flat da imagem de ciência.
        s_f_b_img_fin = f_b_img_fin - sky_bias_flat1 #corrigir o background da imagem de ciência.

        
        hdu = fits.PrimaryHDU() #criando o Header. 
        hdu.data = s_f_b_img_fin
        hdu.header = hdr #header do arquivo corrigido igual ao antigo com algumas flags adicionadas abaixo.

#adicionar os flags dentro do header novo.
        hdu.header['BIAS'] = 'YES'
        hdu.header.comments['BIAS'] = 'correção para bias'
        hdu.header['FLAT'] = 'YES'
        hdu.header.comments['FLAT'] = 'correção para flat'   
        hdu.header['BACKGROUND'] = 'YES'
        hdu.header.comments['BACKGROUND'] = 'correção para background'
        
        orig_file = i.split('.fits') 
        arquivo_final = (x + outfile for x in orig_file)
                
        hdu.writeto(arquivo_final)
  
    






    
     
        




    
















     
    

































    
    
   
    
	
