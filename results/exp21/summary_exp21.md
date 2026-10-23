# Exp21 — Standart bölünmeyi oda sıcaklığına indirme denemesi (SKY130, 300 K)

Tarih: 2026-10-15. Kod: `experiments/exp21_room_temp/rt21.py` (scan, predict, run, predict2, run2, predict3, run3,
analyze1–3, fig). Tahminler koşulardan önce yazıldı: `predictions.json`, `predictions2.json`, `predictions3.json`.
Şekil: `fig_exp21.png`.

## Kısa sonuç

**Olumsuz sonuç; mimarinin sıcaklık sınırı bulundu.** Exp17'deki standart bölünme, indirgenmiş SKY130 modelinde 300 K
değerleriyle oluşmuyor. Toplam 21 cihaz koşusunun hiçbirinde yasanın öngördüğü bölünme yok. Bu bir hata değil, bölünmenin hangi
koşullara dayandığının ölçülmesi. Nedeni üç parçalı:

1. 300 K'de bistabilite zayıf. Pencere 43–68 mV (4 K'de ~770 mV), açık/kapalı akım oranı 5–12 (4 K'de ~10⁴), j_on
   0,08 µA/µm (4 K'de 2,57). Bistabilite yalnız bilinçli eklenmiş bir kuyu sızıntısıyla (g_sub ≈ 1e-11 S/µm) ve kapı
   eşiğin ~0,1 V altındayken var.
2. 300 K'de p-kuyusu donmuyor, iletken. Gerçek kuyuda (R_□ ≈ 3 kΩ/□) gövde bağlanma uzunluğu 2–5 mm, yani 600 µm'lik
   gövde tek düğüm. Cihazın tamamı tek adımda açılıyor ve alan oluşmuyor (R1 ✓).
3. Hücreler bağımsız yapılırsa (ayrı gövdeler) SKY130 uyumsuzluğu ince pencereyi yiyor. 0,5 µm hücre başına σ_VT = 12 mV,
   hücrenin bistabil penceresini ~50 mV kaydırıyor. Sonuç tuz-biber desen. Uyumsuzluk 10 kat azaltılsa bile 1200 hücre
   içindeki en zayıf hücre, bir alan bölünme düşüşüne ulaşmadan çekirdekleniyor. Alanlar bölünerek değil, yeni
   doğumla çoğalıyor (0/13 bölünme).

Varsayımsal çok dirençli bir gövdede (1e8 Ω/□, bağlanma ~10–30 µm) uyumsuzluk ortalanıyor. Tek, düzgün bir alan
çekirdeklenme olmadan büyüyor, ama yasa genişliğinde bölünmüyor; gövde bağlanması merkezin tek başına sönmesini
engelliyor. Bu gövde değeri SKY130'da yok.

4 K bölünmesi şunu istiyor: uyumsuzluğun çok üstünde geniş bir bistabil pencere ve birbirinden bağımsız mikro hücreler
(donmuş kuyu). 300 K'de SKY130'da ikisi de kayboluyor.

## Yöntem

- 4 K modelinin 300 K karşılığı: aynı BSIM4 tabloları (I_d, I_sub, V_T(V_bs)); 4 K düğmeleri 300 K değerine çekildi
  (F_ii = 1, SS 75–95 mV/dek, gövde diyotu V_s = kT/q (26–30 mV), V_on 0,55–0,70 V, g_sub 1e-16…1e-11 S/µm,
  V_ov0 −0,1/0 V).
- Tarama: 288 set, 11 bistabil, 8 Maxwell rejiminde. Yüksek durum V_b ≈ 0,28–0,37 V, yani 4 K'deki "V_b > 0,5 V"
  tanımı geçersiz; yüksek durum kök yapısından tanımlandı.
- Seçilen setler:
  - RT1: SS 85 mV/dek, V_s 26 mV, V_on 0,70 V, g_sub 1e-11, V_ov0 −0,1. İnceltilmiş sabitler: V_down 1,4472,
    V_M 1,4652, V_up 1,4902 V; dV_c 17,96 mV; j_on 0,0789 ve j_low 0,0105 µA/µm.
  - RT2: SS 75 mV/dek. dV_c 25,5 mV; üst pay V_up − V_M = 42,3 mV.
- 300 K'de "kapalı" bölge akım çekiyor (j_low ≈ j_on'un %13'ü). Yasa, kontrast Δj = j_on − j_low ile yazıldı:
  λ_c = √(2 dV_c/(r_D Δj)), w = (I − j_low W)/Δj.

## Sonuçlar ve tahminler

**Tur 1** (`predictions.json`; r_D = 300 Ω/µm, λ = 45–120 µm, R_L = 160 kΩ; 9 koşu):

| tahmin | sonuç |
|---|---|
| P1 filaman oluşur, en az biri kalır | ✓ biçimsel olarak; ama tam uyumsuzlukta 123–145 parçalı tuz-biber desen (en geniş 8–14 µm) |
| P2 toplam genişlik yasası ±%15 | yalnız düşük uyumsuzlukta ✓ (−%11,6); 8 tuz-biber koşusunda ✗ (−%40…−%44) |
| P3 ilk bölünme yasa genişliğinde | ✗ — gerçek bölünme yok; tuz-biberde sayılan "bölünmeler" iz sürme yapaylığı; düşük uyumsuzlukta 7 alan, hepsi çekirdeklenmeyle |
| P4, P5 | değerlendirilemez (gerçek bölünme yok) |

- Kayıtlı tahminlerdeki V_down ve V_up 5 mV'lık kaba ızgaradan geliyordu; dV_c 15,2 yerine 17,96 mV. Tahminler bu
  yüzden de kusurlu, ama sonuç bundan bağımsız olarak olumsuz.

**Tur 2: tasarım penceresi** (`predictions2.json`; sonradan türetilen koşul dV_c < r_v Δj/2 < V_up − V_M; 8 koşu):

| tahmin | sonuç |
|---|---|
| Q1 pencere içi → bölünür | ✗ — 5 koşunun hiçbiri bölünmedi (7–55 çekirdeklenme) |
| Q2 pencere üstü → bölünmeden önce çekirdeklenme | ✓ (ama her durumda böyle) |
| Q3 pencere altı → bölünme yok | ✓ (ama her durumda böyle) |
| Q4 w_split(λ 50)/w_split(λ 30) | değerlendirilemez |

- Neden: 1200 hücrenin uç değeri ~3,2σ. Uyumsuzluk ×0,1'de bile en zayıf hücrenin V_up'ı ~16 mV aşağıda. Bu, RT1'in
  25 mV'luk ve RT2'nin 42 mV'luk payını bölünme düşüşünün altına indiriyor.

**Tur 3: gövde bağlanması** (`predictions3.json`; RT2, tam uyumsuzluk):

| koşul | sonuç |
|---|---|
| R_□ = 3e3 (gerçek kuyu), R1 tahmini | ✓ — cihaz tek adımda tamamen açılıyor (frac_high 0 → 1), en fazla 1 alan |
| R_□ = 1e8 (keşif) | tek alan, çekirdeklenme yok, I ile düzgün büyüyor; yasa genişliğinde (121/202 µm) bölünme yok; ancak 390–534 µm'de asimetrik kopma ve geri birleşme |

## Sonuç ve sınırlar

- Bu modelde ve SKY130'da 300 K bölünmesi yok. Olumsuz sonuç gerekçeleriyle kayıtlı.
- Modelin dışında kalan ve oda sıcaklığında umut veren yol: yüzer gövdeli aygıtın parazitik bipolar kilitlenmesi
  (tek transistör kilitlenmesi, (M−1)β ≥ 1). Bu çok daha güçlü bir bistabilite. 300 K M tablosu 2,6 V'ta bitiyor;
  kilitlenme için V_D ≳ 2,8–3 V gerekir, bu da 1,8 V aygıtın sınırı dışında (5 V'luk kalın oksitli aygıt gerekir).
  Denenmedi.
- Diğer yol: her segmente tasarlanmış bir devre hücresi (mandal) koymak. Bu artık "fiziksel" bölünme olmaz, devre
  bölünmesi olur. Denenmedi.
- Tüm sonuçlar indirgenmiş modele ve 300 K tablolarına bağlı; ölçüm yok.
