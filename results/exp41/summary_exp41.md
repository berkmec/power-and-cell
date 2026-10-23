# Exp41 — Hedefsiz serbest hız laboratuvarı

Tarih: 2026-10-19. Sanal laboratuvar kuruldu; fiziksel laboratuvar ölçümü yapılmadı.

## Kurulum

Bu koşuda `n` veya γ hedefi kullanılmadı. Hız, zaman integrasyonunda serbest durum değişkenidir. γ, ε ve `log10γ` yalnız gösterge/okuma olarak kaydedildi.

**Tezgâh A — sonlu iyon gemisi:** Exp34 EMI fizibil, In diz ve In 5 kV hatları; 10 kg ıslak, 8 kg kuru, 2 kg yakıt, 100 W üst sınırı. Durum değişkenleri kütle, rapidite, uzaklık, gemi öz zamanı ve elektrik enerjisidir. Gaz profilleri tam gaz, %25, rampa ve %50 darbedir. Koşu yalnız yakıt bitince durur.

**Tezgâh B — dış beslemeli gaz tezgâhı:** 8 kg test kütlesi sabit tutulur ve enerji dışarıdan verilir. 0,1g/1g/10g tam gaz ile 1g rampa ve 1g darbe profilleri öz zamanda integre edildi. Bu tezgâh kapalı yakıtlı bir gemi değildir; gereken dış enerji ve anlık güç ayrıca yazılır. Tam-gaz koşularındaki φ=80 sınırı fiziksel hedef değil sayısal güvenlik ufkudur.

Öngörüler koşudan önce `predictions_exp41.json` içinde G1–G6 olarak kaydedildi; altısı da PASS verdi.

## Tezgâh A sonucu: gaz profili son hızı değiştirmiyor

| Hat | Son Δv | Tam gaz | %25 gaz | Rampa | Darbe %50 |
|---|---:|---:|---:|---:|---:|
| EMI fizibil | 0,850 km/s | 0,00615 yıl | 0,02459 yıl | 0,00645 yıl | 0,01229 yıl |
| In diz | 10,001 km/s | 0,711 yıl | 2,844 yıl | 0,747 yıl | 1,422 yıl |
| In 5 kV | 20,455 km/s | 2,963 yıl | 11,850 yıl | 3,111 yıl | 5,777 yıl |

Aynı hatta dört profil arasındaki son-Δv yayılımı yalnız `7,4–7,7×10^-9` kesirdir. Bunun nedeni sürüklenme olmaması ve bütün profillerin aynı 2 kg yakıtı aynı egzoz hızıyla yakmasıdır. Gazı kısmak son hızı düşürmez; yalnız yanış saatini uzatır. En hızlı sonlu hat bile `γ−1=2,33×10^-9` ile relativistik rejimin çok altındadır.

## Tezgâh B sonucu: serbest hız c'ye asimptotik yaklaşır

1g tam gazda olaylar hedeften değil hareket denklemlerinden okundu:

| Okuma | γ | Gemi öz zamanı | Dünya zamanı | Uzaklık | 8 kg kinetik enerjisi |
|---|---:|---:|---:|---:|---:|
| γ=10 | 10 | 2,900 yıl | 9,639 yıl | 8,718 ışık yılı | 6,471×10^18 J |
| In λ=10^-18 m | 11,592 | 3,043 yıl | 11,188 yıl | 10,261 ışık yılı | 7,616×10^18 J |
| In λ=Planck boyu | 7,172×10^17 | 40,499 yıl | 6,948×10^17 yıl | 6,948×10^17 ışık yılı | 5,158×10^35 J |
| γ=10^30 okuması | 10^30 | 67,588 yıl | 9,687×10^29 yıl | 9,687×10^29 ışık yılı | 7,190×10^47 J |

Sabit γ okumasında öz zaman, Dünya zamanı ve uzaklık `1/a` ile ölçeklenir: 0,1g on kat yavaş, 10g on kat hızlıdır; enerji faturası ivmeden bağımsızdır. γ=10^30 anında bile `1−β=5×10^-61`; hız c'ye ulaşmaz.

Serbest koşu son okumaları:

- 0,1g, 85 öz-yıl: γ=3.234, Dünya zamanı 31.325 yıl, enerji `2,324×10^21 J`.
- 1g tam gaz, 77,497 öz-yıl sayısal ufku: γ=`2,770×10^34`, Dünya zamanı `2,684×10^34 yıl`.
- 1g rampa, 85 öz-yıl: γ=`3,820×10^37`.
- 1g darbe %50, 85 öz-yıl: γ=`5,657×10^18`; aynı sürede indiyum göstergesi Planck boyunun 0,127 katına inmiştir.

Dış güç yüksek hızda `m a c` değerine doyar: 0,1g için 2,35 GW, 1g için 23,5 GW ve 10g için 235 GW. Bunlar 100 W gemi kaynağıyla karşılanabilecek değerler değildir.

## Boyut gözlemi

Birlikte hareket eden indiyum test iyonu göstergesi Planck boyunu 1g'de 40,499 öz-yılda geçer. Hareket denklemleri bu noktada düzgün kalır; yeni bir uzamsal boyut, kırılma veya geçiş gözlenmez. Bu yalnız `hν` yarı-klasik okumasının güven alanını terk ettiğini işaretler. Gemi bir bütün olarak tek kuantum değildir; gösterge yalnız test iyonunundur.

## Sınırlar

Sonlu tezgâh Exp34'ün sabit egzoz hızlı, gazla doğrusal kütle akışlı modelini kullanır. Dış tezgâhta enerji kaynağının fiziksel taşıyıcısı, yapısal ivme sınırı, radyasyon, yıldızlararası ortam ve yönlendirme modellenmemiştir. Sonuçlar bir motorun veya laboratuvar düzeneğinin çalıştığını göstermez.

İlk görsel denetimde hız açığı panelindeki kayan-nokta artefaktı rapidite formülüyle giderildi; simülasyon, bağımsız analitik denetim ve figürler yeniden koşuldu.

## Görsel veri ve dosyalar

- `fig_exp41a_finite_lab.png`: sonlu gemide serbest hız, yakıt ve yanış süreleri.
- `fig_exp41b_external_bench.png`: γ, c'ye kalan açık, zaman ayrışması ve dünya çizgisi.
- `fig_exp41c_events_energy.png`: 1g olay çizelgesi, enerji faturası ve gaz profilleri.
- `results/exp41/lab41.json`: iki tezgâhın zaman serileri, olaylar ve laboratuvar yapılandırması.
- `experiments/exp41_free_throttle_lab/free_throttle_lab41.py`: yeniden kullanılabilir sanal laboratuvar simülatörü.
- `experiments/exp41_free_throttle_lab/check41.py`: analitik bağımsız denetim (tümü True).
- `experiments/exp41_free_throttle_lab/fig41.py`: üç figürün üretimi.
