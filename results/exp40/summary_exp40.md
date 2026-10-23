# Exp40 — Denklem 2 (`hν` hattı): ölçek inişi ve boyut denetimi

Tarih: 2026-10-19. Bu çalışma teorik/sayısal bir denetimdir; laboratuvar ölçümü değildir.

## Sorulan soru ve kısa hüküm

Denklem 1 grubunun `h` üzerine kurulu hattı ayrıca çalıştırıldı:

`ε_n = sqrt(1-β²) = m0c²/(hν) = 0,1^n`, `γ=10^n`, `hν=γm0c²`, `λ_h=c/ν=h/(γm0c)`.

Hüküm üç anlamda ayrıdır:

1. **Birim boyutu değişmiyor.** `m0c²` ve `hν` ikisi de `M L² T^-2` boyutundadır; oranları her n için boyutsuzdur.
2. **Uzunluk ölçeği başka rejimlere iniyor.** `λ_h`, elektrozayıf `10^-18 m` çizgisini, Planck enerjisini, kendi Schwarzschild yarıçapını ve Planck uzunluğunu sırayla geçer.
3. **Uzamsal boyut değişimi gözlenmiyor.** Denklemde ek boyut yarıçapı, bağlaşım veya boyut sayısı parametresi yoktur. Planck ölçeğine iniş, ek boyuta giriş kanıtı değildir; yarı-klasik `hν` tanımının standart alanını terk ettiği yerdir.

Öngörüler koşudan önce `predictions_exp40.json` içinde V1–V6 olarak kaydedildi; altısı da PASS verdi. En büyük bağıl sayısal artık `3,33×10^-16` oldu.

## Türlere göre ölçek eşikleri

`λ_C=h/(m0c)` indirgenmemiş Compton dalga boyudur. Eşikler sürekli n değerleridir.

| Tür | λ_C [m] | λ=10^-18 m | E=Planck enerjisi | λ=kendi r_s değeri | λ=Planck boyu |
|---|---:|---:|---:|---:|---:|
| Elektron | 2,426×10^-12 | 6,385 | 22,378 | 22,627 | 23,176 |
| Proton | 1,321×10^-15 | 3,121 | 19,114 | 19,363 | 19,913 |
| İndiyum | 1,159×10^-17 | 1,064 | 17,057 | 17,306 | 17,856 |
| Planck-kütleli test parçacığı | 1,016×10^-34 | −15,993 | 0,000 | 0,249 | 0,798 |

Eşik sırası her noktasal türde aynıdır: önce Planck enerjisi, `0,249 dex` sonra `λ=r_s`, ondan `0,550 dex` sonra `λ=ℓ_P`. Bu farklar `h` ile `ħ=h/(2π)` tanımlarının sabit kaymalarından gelir; yeni bir boyut geçişi değildir.

İndiyum için n=1 değerleri `ν=2,586×10^26 Hz`, `λ_h=1,159×10^-18 m` ve `E=1,070×10^12 eV` olur. Elektron için aynı basamak `5,110 MeV`, proton için `9,383 GeV` verir.

## `hν` boyu ile de Broglie boyu farkı

Aynı γ'ye sahip parçacığın de Broglie dalga boyu `λ_dB=h/(γm0v)` olduğundan `λ_h/λ_dB=β` bulunur. Foton momentumunun parçacık momentumuna oranı `1/β`dır. n=1'de foton-eşdeğeri momentum %0,503 fazladır; n büyüdükçe fark sıfıra gider. Yani Denklem 2'nin pozitif dalı işaret ve birim bakımından tutarlıdır, fakat `hν`yu parçacığın toplam relativistik enerjisine eşitlemek fiziksel bir özdeşleştirme varsayımıdır; sonlu n'de `λ_h`, parçacığın de Broglie boyunun kendisi değildir.

## Makroskopik gemiye uygulama sınırı

8 kg kuru gemi karşıolgusal biçimde tek kuantum sayılsa `λ_h(n=0)=2,763×10^-43 m` olur. Bu, Planck boyunun yalnız `1,71×10^-8` katıdır. Dolayısıyla `hν` hattı geminin tamamına tek-kuantum gibi uygulanamaz; merdiven başlamadan standart yarı-klasik alanın dışındadır. Denklem parçacık/iyon ölçeğinde tanımlanmalıdır.

## Varsayımsal ek boyut katmanı

Ek olarak yarıçapı `R` olan tek kompakt boyut varsayımı denklemin üstüne ayrı bir model olarak konuldu. Geometrik erişim ölçütü `λ_h≤R`, ilk Kaluza–Klein (KK) enerji ölçütü `E≥ħc/R` alındı. Bu iki eşik yarıçap tanımında `log10(2π)=0,798 dex` farklıdır. Örneğin indiyumda `R=10^-18 m` için geometrik eşik n=1,064, toplam-enerji KK eşiği yaklaşık n=0,266 olur. `R=ℓ_P` için değerler 17,856 ve 17,057'dir.

Bu, ek boyutun varlığını göstermez. Denklemin kendisinde `R` olmadığı için hangi n'de "boyuta girileceği" denklemden türetilemez. Böyle bir iddia için bağımsız bir kompaktlaşma modeli ve gözlenebilir dağılım/kuvvet yasası değişimi gerekir.

## Görsel veri

- `fig_exp40a_scale.png`: elektron, proton ve indiyum için dalga boyu, frekans ve enerji merdivenleri; `λ_h`–de Broglie farkı.
- `fig_exp40b_dimension_audit.png`: kütle–γ rejim haritası, Planck boyuna ve Schwarzschild yarıçapına iniş, tür bazında eşik çubukları.
- `fig_exp40c_extra_dimension.png`: varsayımsal `R` için geometrik/KK eşikleri, KK aralığı, mod kapasitesi ve indiyum `R–n` erişim haritası.

Üç figür de üretildikten sonra görsel olarak açıldı; çakışan açıklama ve kontur etiketi biçimleri düzeltilip yeniden kontrol edildi.

## Dosyalar

- `results/exp40/predictions_exp40.json` — koşudan önce kayıtlı V1–V6
- `results/exp40/h40.json` — beş kütle için 0–30 sürekli merdiven, eşikler ve KK katmanı
- `results/exp40/fig_exp40a_scale.png`
- `results/exp40/fig_exp40b_dimension_audit.png`
- `results/exp40/fig_exp40c_extra_dimension.png`
- `experiments/exp40_hnu_dimension/h40.py` — simülasyon
- `experiments/exp40_hnu_dimension/check40.py` — CODATA ve satır-kökü bağımsız denetimi (tümü True)
- `experiments/exp40_hnu_dimension/fig40.py` — üç figürün üretimi
