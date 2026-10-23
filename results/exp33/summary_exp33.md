# Exp33 — SKY130 hücresinin 4 K'de ölçümü: mevcut veriden çıkarılabilen her şey + laboratuvar kaydı

Tarih: 2026-10-16. Kod: `experiments/exp33_cell_4k/`
- `cell33.py`: C33-1…C33-6 (kayıtlı ölçütler);
- `cell33_sens.py`: kayıt-sonrası keşif (basamak tanısı + gövde-etkisi duyarlılığı) — hüküm değildir.
- `fig33.py`: `fig_exp33.png` istatistik figürü (4500×1625 px).

Öngörüler hesaptan önce yazıldı: `predictions_exp33.json` (hedef gözlenebilirler hesaplanmadan önce; dayanak açıkça kör-olmayan olarak bildirildi).

Sonuçlar:
- `cell33.json` (envanter + C33-1…5 + Part B laboratuvar tablosu);
- `cell33_sensitivity.json` (keşif);
- `fig_exp33.png` — gövde etkisi, kink yükseltimi ve basamak denetimi panelleri.

**Fiziksel laboratuvar ölçümü yapılmadı.** Bu bilgisayarda kriyostat, balastlı geniş test yongası ve ölçü aleti yoktur; eldeki tek 4 K ölçümü NIST mds2-2997 ham verisidir (en geniş aygıt 7 µm < w_s = 20–54 µm, balastsız, yalnız ileri tarama). Bu deney iki şey yapar: (A) ham 4 K veriden Exp32'nin çıkarmadığı hücre girdilerini ölçer; (B) geniş/balastlı laboratuvar koşusunu sayı ve protokolle kaydeder (C33-6/7).

## Kısa sonuç

- **C33-1 ✗ (gövde etkisi).** VB yalnız ters yönde (0, −0,75, −1,5 V). Standart aygıtlarda ΔVT(−1,5 V) = +0,147…0,237 V, |dVT/dVB| = 0,098…0,158 (medyan 0,105). Exp12'den beri `dvt_of` içinde kullanılan BSIM4 300 K çapası aynı büyüklükte 0,316 V ve ortalama eğim ~0,32 verir: ölçülen medyan oran 0,50, eğim oranı ~0,33. Yani **4 K gövde etkisi model varsayımının ~2–3 kat altında**; işaret 10/10 doğru. Uyarı: ölçülen ters yöndür, model iç gövde yükselişinde (ileri yön) aynı büyüklüğü kullanıyordu; vekil karşılaştırma öngörüde açık yazıldı. Sonuç: hücre sabitlerindeki iki varsayımdan biri ölçüldü ve **tutmadı**; ileri yön hâlâ ölçülmedi.
- **C33-2 ölçülemez (gövde diyotu).** 4 K IDVG dosyalarında IB sütunu yoktur (başlıklar ID/IG/VG + türetilmiş sütunlar) ve VB > 0 bloğu yoktur. Von = 1,05 V, Vs = 5 mV varsayımı aynen kalır.
- **C33-3 ✓ (kink yükseltimi).** A = ΔID_fazla/|IB| @1,8 V, VG = 0,9 V: medyan **1544** (n = 8, aralık 475–6333, hepsi > 30). Fazla drenaj akımı gövde akımının kendisi değil, ~1,5·10³ katı eşik-kayması yükseltimidir — Exp12 donmuş-gövde mekanizması 10 aygıtın tamamında nicel olarak doğrulandı.
- **C33-4 kayıtlı biçimiyle ✗, 4 K alt iddiası ✓ (hücre basamağı yokluğu).** 4 K'de 44 IDVD taraması: 0,5 dex üstü bitişik-nokta basamağı **0**, en büyük 0,386 dex (lvt W0,42/L0,15, VG = 0,6 V, ID 13,7→33,2 nA @VD = 1,625 V: eşik-altı kink sürekliliği, µA–mA kilitlenmesi değil). Hüküm yine de fail yazılır çünkü kayıtlı ölçüt 300 K kontrolünü de kapsıyordu ve orada tek 0,641 dex basamak var (W7/L8, VG = 0, ID 4,52·10⁻¹⁰→1,98·10⁻⁹ A @VD = 0,9 V — kapalı-durum sızıntı tabanı, alet artefaktı). Fizik iddiası — balastsız dar aygıtta kuantumlu hücre basamağı yok — 4 K'de tuttu.
- **C33-5 doğrulandı (histerezis bu veriden ölçülemez).** 24 dosyanın bütün VD/VG taramaları tek yönlü monotondur; yukarı/aşağı çifti yoktur.
- **Part B laboratuvar tablosu (model, Exp32 ölçülmüş girdiler; Exp32 ızgarasıyla birebir, sapma 0).** V_M = 0,701 V, V_down = 0,639 V, V_up = 2,323 V, dV_c = 0,062 V, j_on = 64,06 µA/µm. Tasarım noktası (r_D = 10, λ = 60 µm): w_s = 31,6 µm, I_c = 1,52 mA, basamaklar 1,52/3,04/6,08/10,64 mA, 4 hücre için W ≥ 127 µm (ızgara geneli 82–217 µm). **V_up = 2,32 V, nfet_01v8'in 1,8 V anmasının üstündedir**: anma içinde kendiliğinden doğum beklenmez, tohum şarttır (Exp18/Exp32 ile tutarlı); histerezis, tohumla doğum + aşağı taramada V_down = 0,639 V ile sınanır.
- **C33-7 güç.** Tasarım noktasında hücre başına V_M'de 1,07 mW, 7 hücre 7,46 mW; ızgara geneli 7 hücre 4,84–12,76 mW @V_M, 1,8 V'da hücre başına 2,74 mW'a kadar. 4 K sahnesinde darbeli tarama / düşük görev çevrimi protokol şartıdır.

## Kayıt-sonrası keşif (hüküm değil): gövde-etkisi duyarlılığı

`cell33_sens.py`, nist32 `measured_transfer` kod yolunda `dvt_of` ölçeğini s ile çarpar (s = 1 kayıtlı sabitleri birebir verir: V_M 0,701012):

| s | V_M | V_down | V_up | j_on | w_s (λ60, r_D10) | I_c |
|---|---|---|---|---|---|---|
| 1,00 (BSIM çapası) | 0,701 V | 0,639 V | 2,323 V | 64,06 µA/µm | 31,6 µm | 1,52 mA |
| 0,50 (ölçülen |ΔVT| oranı medyanı) | 0,934 V | 0,857 V | 2,384 V | 6,53 µA/µm | 198,9 µm | 0,97 mA |
| 0,33 (ölçülen eğim oranı) | 1,149 V | 1,053 V | 2,419 V | 0,85 µA/µm | ∞ (λ ≤ λ_c, bölünme yok) | — |

Sabitler gövde etkisine üstel duyarlıdır çünkü çalışma noktası ölçülen transfer eğrisinin dik omzundadır. Vekil uyarısı geçerlidir: ters ölçüm ileri ölçek sayılmıştır; ileri yön karekök yasasıyla daha dik olabilir. Bu yüzden yeni sabit yayımlanmaz; ancak **ileri-yön gövde etkisi artık laboratuvar öncesi en kritik ölçülmemiş girdidir** ve mevcut ETest-tipi dar yapılarda (yeni üretim gerekmeden) ölçülebilir.

## Laboratuvar protokolü (kayıtlı, C33-6/7)

- **Yapı:** nfet_01v8, W ≥ 150–220 µm (tasarım noktası ≥127 µm), dilimli drenaj + dağıtık poli balast (λ = 40/60/120 µm varyantları, Exp13/17 geometrisi, r_D = 5–20), ayrı gövde kontağı, tohum için dilimli kapı.
- **Ortam/alet:** 4 K kapalı-devre prob istasyonu (NIST düzeneği sınıfı), 4156B/B1500 sınıfı SMU, IB kaydı zorunlu, yıkıcı kilitlenmeye karşı akım uyumu (compliance).
- **Sıra:**
  0. *Önce ileri gövde (C33-1/2 kapanışı):* dar standart aygıtta IDVG, VB = 0…+0,8 V (IB uyumu ile) → ileri dVT(Vb) eğrisi + diyot Von.
  1. Karakterizasyon: IDVD/IDVG çift yön (yukarı+aşağı), VD ≤ 1,8 V (V_up anma dışında; istenirse ayrı kalın-oksit kol).
  2. Tohumlu doğum: +30 mV yerel kapı dilimi, 10 denemede ≥8/10 yer seçimi (dilim akımlarıyla).
  3. Basamaklar: ani V_DD adımıyla k = 1, 2, 4 platoları k·I_c ±%25; zorla düşüşte 3→7.
  4. Ölüm/aktarım: bir hücre öldürülünce toplam akım değişimi ≤%3.
  5. Sıcaklık: 4/20/30/40 K — kink ve kilitlenme 30–40 K arasında kaybolmalı.
  6. Güç disiplini: toplam >10 mW ise darbeli tarama, sahne sıcaklığı kaydı.
- **Riskler:** yıkıcı latch-up (balast+uyum), öz-ısınma (mW mertebesi @4 K sahte V_down kayması üretir), ESD.

## Son duruma etkisi

Exp0–27 sayımı değişmez; Exp33, Exp28–32 gibi sayım dışı kapanış/ölçüm deneyidir. Hücre sabitleri (V_M = 0,70 V, j_on = 64 µA/µm, I_c = 1,0–2,6 mA) ileri-gövde ölçümüne kadar "BSIM gövde çapası varsayımıyla" etiketi taşır; vekil duyarlılık aralığı V_M 0,70–1,15 V ve j_on 0,85–64 µA/µm'dir. Açık kalan tek fiziksel iş hâlâ aynıdır: geniş/balastlı SKY130'da kilitlenme + histerezisin kendisi — fakat artık protokolü, sayısal hedefleri ve öncelikli ön-ölçümü (ileri gövde) kayıtlıdır.

## Kaynaklar

- NIST mds2-2997 ham .mdm (4 K/300 K, `results/exp32/nist_sky130/`) — Exp32'deki kopya, dokunulmadı.
- Exp12 `reduced_bistability.json` BSIM4 gövde çapaları; Exp17 `law.py`; Exp32 `nist32.py` tanımları ve `nist_sky130.json` sabitleri.
