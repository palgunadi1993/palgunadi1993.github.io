---
layout: post
title:  "'Bootstrap' uncertainty analyis for hypoDD (Bahasa)"
date:   2019-12-21 00:00:00
categories: software manual
---

Setelah melakukan relokasi hiposenter dengan menggunakan teknik *double-difference* akankah lebih baik jika mencantumkan error 3D (x , y, z). Secara umum error dari lokasi ini direpresentasikan dengan 3D elips. Namun, *uncertainty* spasial dari hiposenter tidak dapat diperoleh secara langsung dari output hypoDD. Untuk mendapatkan spasial error hanya dapat dihitung secara statistik. Terdapat berbagai jenis analisis statistik yang dapat dilakukan, salah satunya adalah metode 'bootstrap'. Pada tutorial kali ini, saya hanya akan membahas metode statistik bootstrap untuk mengidenditifikasi spatial error dari hasil relokasi mengunggunakan teknik *double-difference*. 

Sebelumnya, jika ingin menggunakan metode ini pada relokasi hypoDD, jangan lupa untuk mensitasi "Supendi, P., Nugraha, A.D., Widiyantoro, S. et al. Geosci. Lett. (2019) 6: 18. https://doi.org/10.1186/s40562-019-0148-9" ([link][paper]). Terima kasih.

Teori dasar tentang metode ini telah dijelaskan pada paper tersebut diatas. Disini, saya akan lebih membahas secara teknis. Sebelum menjalankan script berikut, pastikan dahulu anda telah memiliki program hypoDD yang dapat di download secara gratis di website USGS. 

Step-by-step yang harus dilakukan adalah:
1. Set parameters dari hypoDD dan run relokasi.
2. Setelah anda yakin dengan semua parameters tersebut dan hasil dari relokasi memuaskan, copylah script berikut dan letakkan didalam folder hypoDD_folder/src/ 

{% highlight python %}
import os
import numpy as np
import time

# manipulating earthquake data
fname = 'ph2dt/PaluEQ.pha'
data = open(fname)
with open(fname) as f:
    content = f.read().splitlines()

firstLine = data.readlines()
for i in range(len(firstLine)):
	firstLine[i]=firstLine[i].split()
data.close()

#standard deviation of pick-time in seconds
std = 0.1

# number of realization
N = 1000

for nn in range(N):
    #write DATA
    hypodd = open('ph2dt/stats/PaluEQ'+str(nn)+'.pha', 'w+')
    # hypodd.write(content[0]+'\n')
    j = 1
    while j < len(content)-1:
        hypodd.write(content[j - 1]+'\n')
        i = 0
        while firstLine[j][0] != '#':
            i += 1
            j += 1
            if j == len(content):
                break
        j += 1

        # Generating random number with normal distribution
        rand = np.random.normal(loc = 0, scale = std, size = i)

        k = j-i-1
        count = np.arange(k, j - 1, 1)
        r = 0
        for m in count:
            tp = np.round(float(firstLine[m][1]) - rand[r], decimals=2)
            hypodd.write('%5s'%firstLine[m][0]+'%12s'%str(tp)+'%8s'%firstLine[m][2]+'%4s'%firstLine[m][3]+'\n')
            r += 1

    hypodd.close()

    time.sleep(2)
    # modify and run ph2dt
    print('===================================')
    print('===================================')
    print('=========Now running ph2dt=========')
    print('===================================')
    print('===================================')
    with open('ph2dt/stats/ph2dt.inp', 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    # now change the line
    data[4] = 'stats/PaluEQ'+str(nn)+'.pha\n'

    # and write everything back
    with open('ph2dt/stats/ph2dt.inp', 'w') as file:
        file.writelines( data )
    cmd = 'cd ph2dt && ./ph2dt stats/ph2dt.inp'
    os.system(cmd)

    time.sleep(3)
    print('===================================')
    print('===================================')
    print('=========Now running hypoDD========')
    print('===================================')
    print('===================================')

    with open('hypoDD/bootstrap/hypoDD.inp', 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    # now change the line
    data[18] = 'hypoDD/bootstrap/hypoDD'+str(nn)+'.reloc\n'

    # and write everything back
    with open('hypoDD/bootstrap/hypoDD.inp', 'w') as file:
        file.writelines( data )


    cmd = 'hypoDD/hypoDD hypoDD/bootstrap/hypoDD.inp'
    os.system(cmd)
    time.sleep(1)
{% endhighlight %}
 
3. Pastikan anda memilih standar deviasi dari picking dan jumlah realisasi yang sesuai. Besarnya jumlah realisasi akan sangat berpengaruh dengan lamanya analisis yang akan dilakukan. Be wise :)

4. Untuk plotting dapat menggunakan script berikut (letakkan pada folder hypoDD_folder/src/):
{% highlight python %}
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.basemap import Basemap

N = 1000
M = 567
hypo = np.loadtxt('hypoDD/hypoDD.reloc', usecols=(0,1,2,3))
id = np.array([int(l)-1 for l in hypo[:,0]])
hx = hypo[:,2]
hy = hypo[:,1]
hz = hypo[:,3]

# initialize matrix
X = np.zeros((M,N))
Y = np.zeros((M,N))
Z = np.zeros((M,N))

for i in range(N):
    data = np.loadtxt('hypoDD/bootstrap/hypoDD'+str(i)+'.reloc', usecols=(0,1,2,3))
    idx = np.array([int(l)-1 for l in data[:,0]])
    X[idx,i] = data[:,2]
    Y[idx,i] = data[:,1]
    Z[idx,i] = data[:,3]

X[X==0]=9999
Y[Y==0]=9999
Z[Z==0]=9999

## MapView
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
relative_error = open('relative_error.txt', 'w+')
relative_error.write('x_err'+'\t'+'y_err'+'\t'+'z_err'+'\n')
ll = 0
for mm in id:
    x = X[mm,:]
    y = Y[mm,:]
    z = Z[mm,:]

    x = np.delete(x, np.where(x==9999))
    y = np.delete(y, np.where(y==9999))
    z = np.delete(z, np.where(z==9999))

    cov = np.cov(x,y)
    lambda_, v = np.linalg.eig(cov)
    lambda2 = np.sqrt(lambda_)

    cov = np.cov(x,z)
    lambda_, v = np.linalg.eig(cov)
    lambdaz1 = np.sqrt(lambda_)

    cov = np.cov(y,z)
    lambda_, v = np.linalg.eig(cov)
    lambdaz2 = np.sqrt(lambda_)

    tempx = lambda2[0] * 111.11
    tempy = lambda2[1] * 111.11
    tempz = (lambdaz1[1] + lambdaz2[1])*0.5

    relative_error.write(str(tempx)+'\t'+str(tempy)+'\t'+str(tempz)+'\n')

    ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                      width=lambda2[0]*2*2, height=lambda2[1]*2*2,
                      angle=-np.rad2deg(np.arccos(v[0, 0])))
    ell.set_facecolor('none')
    ell.set_edgecolor('k')
    ell.set_linewidth(0.5)
    ax.add_artist(ell)
    ax.scatter(np.mean(x), np.mean(y), 0.1, 'k')
    # ax.scatter(x,y,2, marker='+', color='b')
    ll += 1
relative_error.close()
# ax.set_xlabel('Longitude [$^o$]')
# ax.set_ylabel('Latitude [$^o$]')
ax.set_xlim([119, 121.5])
ax.set_ylim(-3, 1)
m = Basemap(llcrnrlat=-3,urcrnrlat=1,\
            llcrnrlon=119,urcrnrlon=121.5,lat_ts=1,resolution='h', ax=ax)
m.drawcoastlines(linewidth=0.4)
m.fillcontinents(color='lightgray',lake_color='white')
m.drawmapboundary(fill_color='white')
m.drawparallels(np.arange(-3,1,1),labels=[1,0,0,0], linewidth=0.0)
m.drawmeridians(np.arange(119,121.5,1),labels=[0,0,0,1], linewidth=0.0)
fig.savefig('MapView_bootstrap1.png', fmt='png', dpi=500, bbox_inches='tight')
plt.show()
{% endhighlight %}

5. Jika terdapat error pada file python diatas, berarti anda belum meng-install packages yang diperlukan. Untuk mendapatkan packages tersebut, silahkan mencari informasinya di google ;).

Sekian tutorial kali ini, semoga bermanfaat. 

> Jika ada pertanyaan silahkan ajukan ke saya -> hendrawan.palgunadi@gmail.com

Terima kasih.


[paper]: https://www.researchgate.net/publication/337947322_Hypocenter_relocation_of_the_aftershocks_of_the_Mw_75_Palu_earthquake_September_28_2018_and_swarm_earthquakes_of_Mamasa_Sulawesi_Indonesia_using_the_BMKG_network_data
