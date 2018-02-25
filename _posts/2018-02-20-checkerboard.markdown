---
layout: post
title:  "*Checkerboard* SIMULPS12 (Bahasa Indonesia)"
date:   2017-02-25 00:00:00
categories: software manual
---

Pada postingan kali ini saya akan menguraikan cara untuk melakukan tes resolusi dari tomogram yang telah kita lakukan sebelumnya. Namun sebaiknya *checkerboard resolution test* ini dilakukan sebelum menentukan parameter yang akan digunakan pada inversi tomografi. Secara sederhana tes resolusi ini penting untuk mengetahui daerah dengan hasil inversi yang baik. Sebagai contoh adalah daerah yang dilalui oleh sinar seismic (*raypath*), semakin banyak *raypath* yang melalui daerah tersebut akan menghasilkan hasil inversi yang lebih baik dari daerah yang sedikit dilewati *raypath*. Sebagai ilustrasi dapat dilihat pada publikasi konferensi saya [disini][paper]. Namun untuk mengetahui lebih jauh tentang tes resolusi dapat membaca buku karangan (alm) Prof. Albert Tarantola [disini][tarantola], karena untuk tes resolusi dapat menggunakan berbagai macam cara seperti: *matrix identity analysis*, *derivative weight sum*, *diagonal matrix element*, dan lain-lain. Untuk lebih jelasnya dapat membaca pada jurnal-jurnal internasional.

Disini saya hanya akan membahas tes resolusi yang biasa digunakan yaitu *checkerboard resolution test* (sebut saja CKB) dengan menggunakan simulps12. Bagi yang belum memiliki program simulps12.f dapat mengunduhnya pada postingan saya yang sebelumnya ([tautan][simulps]). Ide dasar CKB dengan menggunakan simulps12.f adalah dengan membuat waktu tempuh sintetik kecepatan awal papan catur, selanjutnya dengan menggunakan waktu tempuh papan catur tersebut, kita gunakan sebagai waktu tempuh layaknya data observasi. Hal ini berguna untuk mengetahui daerah yang menghasilkan hasil inversi yang mirip dengan kecepatan awal papan catur tadi. 

### Model kecepatan CKB
Sebelum membuat waktu tempuh sintetik, terlebih dahulu membuat model kecepatan papan catur yang disimpan dengan nama file `velocity.dat`. Berikut adalah *script* yang bisa digunakan untuk membuatnya. 

{% highlight python %}
# -*- coding: utf-8 -*-
"""
@author: palgunadi
"""

nx = 15
ny = 19
nz = 14
nzz = 2*nz

## Kecepatan diikuti oleh nilai vp/vs
vel = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 3.1, 4.2,  4.3, 5.4, 5.5, 5.6,  5.7, 5.8,
       1.73, 1.73,  1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73] 

data = '/path/to/file/velocity.dat'

fileinput = open(data, 'w')
fileinput.write(' 1.0 '+str(nx)+' '+str(ny)+' '+str(nz)+'\n')


## Grid dalam km
#=====Yang diedit untuk grid===================================================
fileinput.write('-200.0 -25.0 -15.0  -9.0  -6.0  -4.0  -2.0   0.0   2.0   4.0   6.0   9.0  15.0  21.0 200.0\n')
fileinput.write('-200.0 -25.0 -20.0 -15.0 -10.0  -8.0  -6.0  -4.0  -2.0   0.0   2.0   4.0   6.0   8.0  10.0  12.0  15.0  25.0 200.0\n')
fileinput.write('-200.0  -2.0  -1.5  -1.0  -0.5   0.0   0.5   1.0   1.5   2.0   3.0   5.0  15.0 200.0\n')
fileinput.write(' 0  0  0  \n')
#==============================================================================

i = 0
j = 0
k = 0
for i in range(nzz):
    if (i+1)%2 == 0:
        for j in range(ny):
            if (j+1)%2 == 0:
                for k in range(nx):
                    if (k+1)%2 == 0:
                        v = vel[i] + (0.1*vel[i]) ## 0.1 menunjukkan nilai 10% dari kecepatan awal
                        fileinput.write(('%5.2f')%v)
                    else:
                        v = vel[i] - (0.1*vel[i])
                        fileinput.write(('%5.2f')%v)
                fileinput.write('\n')
            else:
                for k in range(nx):
                    if (k+1)%2 == 0:
                        v = vel[i] - (0.1*vel[i])
                        fileinput.write(('%5.2f')%v)
                    else:
                        v = vel[i] + (0.1*vel[i])
                        fileinput.write(('%5.2f')%v)
                fileinput.write('\n')
    else:
        for j in range(ny):
            if (j+1)%2 == 0:
                for k in range(nx):
                    if (k+1)%2 == 0:
                        v = vel[i] - (0.1*vel[i])
                        fileinput.write(('%5.2f')%v)
                    else:
                        v = vel[i] + (0.1*vel[i])
                        fileinput.write(('%5.2f')%v)
                fileinput.write('\n')
            else:
                for k in range(nx):
                    if (k+1)%2 == 0:
                        v = vel[i] + (0.1*vel[i])
                        fileinput.write(('%5.2f')%v)
                    else:
                        v = vel[i] - (0.1*vel[i])
                        fileinput.write(('%5.2f')%v)
                fileinput.write('\n')

fileinput.close()

{% endhighlight %}

### File stasiun
Untuk file stasiun, format file dan koordinat sama seperti yang telah saya tulis pada post sebelumnya ([link][simulps]). Tidak perlu ada yang dirubah dari file ini.

### File control
Pembuatan file control hampir sama dengan apa yang telah saya tulis pada postingan sebelumnya, yang membedakan hanya variabel nitloc dan nitmax __diubah__ menjadi __-1__ seperti yang ditunjukkan pada gambar dibawah ini:

![control.dat ckb]({{ "/assets/images/ckb_control_dat.png" | absolute_url }})

### Jalankan ckb untuk menghasilkan waktu tempuh sintetik.
Sebelum mulai menjalankan program (_compiling_), terlebih dahulu pastikan semua file input (`control.dat`, `velocity.dat`, `earthq.dat`, `stations.dat`, `simulps12.f`, `simulps_common.inc`, dan `velfile.33tq1c`) berada pada 1 folder yang sama.

Compile semua file menjadi 1 (dapat dioperasikan pada ubuntu 16.04 dengan gfortran, gcc 5.4.0).
{%highlight bash%}
~/path/to/file/ gfortran -std=legacy simulps12.f -o out.file
{% endhighlight %}
Setelah semua file tersebut ter-compile maka jalankanlah file eksekusi `out.file`
{% highlight bash %}
~/path/to/file/ ./out.file
{% endhighlight %}

## Menjalankan sintetik model
Setelah pemodelan kedepan dijalankan, akan dihasilkan file dengan nama `fort.24` yang merupakan file waktu tempuh sintetik. Selanjutnya file ini akan digunakan untuk melakukan inversi tomografi untuk mengetahui apakah hasil inversi tersebut baik (kembali ke papan catur) atau tidak terjadi update kecepatan. Namun, file ini tidak bisa digunakan secara langsung sebagai waktu tempuh sintetik, melainkan harus diedit kembali agar tidak terjadi *error*. Berikut adalah *script* yang telah saya buat untuk merubah sepenuhnya file `fort.24` menjadi file `earthq.dat`.

{% highlight python %}

## Membaca file fort.24
fileinput = '/path/to/file/named/fort.24'

file = open(fileinput, 'r')
baris = file.readlines()
file.close()

## File earthq.dat
fileoutput = '/path/to/file/earthq.dat' 
earth = open(fileoutput,'w')

for i in range(len(baris)):
    a = baris[i][0:4]
    b = baris[i][14:18]
    c = baris[i][28:32]
    d = baris[i][42:46]
    e = baris[i][56:60]
    f = baris[i][70:74]
    if i == 0:
        earth.write(baris[i][:50]+'\n')
    elif baris[i][0] == '1' or baris[i][0] == '0':
        earth.write(baris[i])
    if len(baris[i]) == 85:
        if a == b:
            earth.write(baris[i][:5]+'P')
        else:
            earth.write(baris[i][:5]+'P')
        if b == c:
            earth.write(baris[i][6:19]+'P')
        else:
            earth.write(baris[i][6:19]+'S')
        if c == d:
            earth.write(baris[i][20:33]+'P')
        else:
            earth.write(baris[i][20:33]+'S')
        if d == e:
            earth.write(baris[i][34:47]+'P')
        else:
            earth.write(baris[i][34:47]+'S')
        if e == f:
            earth.write(baris[i][48:61]+'P'+baris[i][62:75]+'S'+baris[i][76:])
        else:
            earth.write(baris[i][48:61]+'S'+baris[i][62:75]+'P'+baris[i][76:])
    if len(baris[i]) == 57:
        if a == b:
            earth.write(baris[i][:5]+'P')
        else:
            earth.write(baris[i][:5]+'P')
        if b == c:
            earth.write(baris[i][6:19]+'P')
        else:
            earth.write(baris[i][6:19]+'S')
        if c == d:
            earth.write(baris[i][20:33]+'P')
        else:
            earth.write(baris[i][20:33]+'S')
        if d == e:
            earth.write(baris[i][34:47]+'P'+baris[i][48:])
        else:
            earth.write(baris[i][34:47]+'S'+baris[i][48:])
    if len(baris[i]) == 29:
        if a == b:
            earth.write(baris[i][:5]+'P')
        else:
            earth.write(baris[i][:5]+'P')
        if b == c:
            earth.write(baris[i][6:19]+'P'+baris[i][20:])
        else:
            earth.write(baris[i][6:19]+'S'+baris[i][20:])
earth.close()
{% endhighlight %}

Setelah file `earthq.dat` terbentuk, gunakan model kecepatan sebenarnya (bukan model papan catur) sebagai file `velocity.dat`, untuk file `stations.dat` tetap menggunakan file yang sebelumnya. Untuk menjalankan program, gunakan cara yang telah saya sebutkan pada subbab sebelumnya. Akhirnya, akan terbentuk file `fort.23` yang siap untuk di plot.
Selamat mencoba.


> Semoga apa yang saya tulis ini berguna untuk kalangan riset. Jika ada pertanyaan silahkan untuk menghubungi hendrawan.palgunadi@gmail.com


[tarantola]: http://www.ipgp.fr/~tarantola/Files/Professional/Books/InverseProblemTheory.pdf
[paper]: https://www.researchgate.net/publication/316443241_Steam_and_Brine_Zone_Prediction_around_Geothermal_Reservoir_Derived_from_Delay_Time_Seismic_Tomography_and_Anisotropy_Case_Study_PR_Geothermal_Field
[simulps]: https://palgunadi1993.github.io/software/manual/2017/02/20/first.html