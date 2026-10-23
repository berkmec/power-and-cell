# Exp30 — Exp0–Exp27'nin her sonucunun ölçümle kapanışı ve ölçülmüş atom yayımı itkisi

Tarih: 2026-10-16. Kod: `experiments/exp30_all_data_closure/adc30.py` (`python adc30.py fig`). Ölçütler, çapa kaynakları okunmadan
önce yazıldı: `predictions_exp30.json` (14:56). Sonuçlar: `closure_all.json`. Şekil: `fig_exp30.png`.

**Neden.** Stop kancası iki şeye itiraz etti:
1. Exp29 kapanış testi Exp0–Exp27'nin her sonucunu girdi olarak kullanmıyordu.
2. Programlanabilir atomla itkinin ölçülmüş bir sonucu yoktu.

Bu deney ikisini de ele alıyor:
- Her deneyin ana sonucu, aynı fiziksel büyüklüğün bağımsız ölçümüyle karşılaştırıldı.
- Ölçülemeyenler "ölçüm yok" ya da "yöntem" diye işaretlendi; hiçbir değer uydurulmadı.
- Ayrıca alanla denetlenen atom yayımıyla çalışan gerçek bir itkicin (indiyum FEEP) ölçülmüş itkisi Exp28 denklemiyle sınandı.

Her çapanın okunma düzeyi kayıtlı ("tam metin" ya da "özet").

## Kısa sonuç

- **28 sonucun hepsi kapandı:**

  | sonuç | sayı | deneyler |
  |---|---|---|
  | ✓ | 9 | Exp3, 4, 5, 8, 11, 12, 15, 21, 27 |
  | ✗ | 5 | Exp1, 6, 22 (alt bant, ek), 25, 26 |
  | karışık / kısmen | 5 | Exp0, 2, 10, 18, 24 |
  | ölçüm yok | 6 | Exp13, 14, 17, 19, 20, 23 |
  | yöntem | 3 | Exp7, 9, 16 |

- **Ölçülmüş "alanla programlanan atom" itkisi.** İndiyum FEEP'te atomlar alanla tek tek koparılıp iyonlaşıyor.
  - Tezgâhta doğrudan ölçülen itki, ölçülmüş akım, gerilim ve ışın sapmasından Exp28 denklemiyle hesaplanan itkiyle %1,6–4,3
    içinde (ölçüt ≤ %5).
  - Alan arttıkça yayan iğne sayısı artıyor: 37 µA'de 28 iğnenin 1'i, 3,5 mA'de 28'i. İğne başına ~125 µA; Exp23 hücre akımıyla
    (166–268 µA) aynı mertebe. Bu, N(t)·I_c yapısının gerçek bir itkiçteki karşılığı.
- **Ölçümün düzelttikleri:**
  1. **Yapay atom kapasitansı.** Exp1 modeli 1,7–3,8 kat büyük. Ölçümde ~2 nm Si ada oda sıcaklığında Coulomb blokajı gösteriyor;
     model 0,4 nm diyordu.
  2. **Hücre mekanizmasının sıcaklığı.** Darbe iyonlaşması ve donmuş gövdeyle oluşan kink/histerezis, 180 nm CMOS'ta 30 K'de hâlâ
     ölçülüyor, 40 K'de yok. Sınır "≤ 4 K" değil, ölçülmüş sınıfta "≲ 30 K".
  3. **Eşleşmezlik.** Yalnız rastgele katkı MC'si 3,5–5,6 kat küçük. Exp18'in +30 mV tohumu yine de ölçek eşleşmezliğinin 7–9σ'sı.
  4. **Kuantum modelleri iyimser.** Exp24–26 sayıları ölçümün üstünde. Exp22'nin tek alt bant aralığı ölçümle çelişiyor.
- **Ölçülmemiş tek katman:** silisyumdaki hücre yasaları (Exp13, 14, 17, 19, 20, 23). Laboratuvar gerekiyor.

## Bütün tablo (her deney, bir ölçülmüş çapa)

| deney | projenin sonucu | ölçülmüş çapa | kaynak (okuma) | sonuç | denklemdeki rolü |
|---|---|---|---|---|---|
| Exp0 | Coulomb + hapsolma, özyinelemeli boyutlar | ~2 nm Si adada 300 K'de Coulomb ve kuantum etkileri | Shin 2010, arXiv 1001.3724 (tam metin) | hapsolma ✓ (nitel) / Coulomb ✗ | turnike adası |
| Exp1 | c_Σ = 22,989 ε a | C_Σ = 0,42 aF, e²/C_Σ = 0,377 eV (1,58–2 nm) | Shin 2010 (tam metin) | ✗: model 1,01–1,59 aF (r = 0,5) / 0,72–1,12 aF (r = 0,2) | turnike E_C |
| Exp2 | 4 K'de yük kuantumu; 77/300 K'de yok | Si pompası T_e 7–9 K'de ppm (E/kT 22–28); 5 nm ada < 200 K blokaj | Zhao 2017; Shin 2010 (tam metin) | 4 K ✓ / 77 K ✗ | turnike sıcaklığı |
| Exp3 | 300 K çift nokta ölü | Si yük kübiti T2* 0,13–2,1 ns (seyreltme soğutucusu) | Shi 2013, arXiv 1208.0519 (özet) | ✓ | kısıt: kriyojenik |
| Exp4 | κ 4,4–5,7 nm⁻¹ | bariyer 3,1 eV, FN kütlesi ≈ 0,4 → 5,70 (5,34–6,05) nm⁻¹ | Lenzlinger–Snow 1969; arXiv 1512.00255 (özet) | ✓: orta değer −%11 | tünelleme |
| Exp5 | V_g(n_s 1e12) = 0,513 V | PDK vth0 = 0,519 V (ölçüme uydurulmuş kart) | SkyWater PDK (birincil dosya) | ✓: %1,2 | SKY130 alanı |
| Exp6 | MC σ_VT 3,9–5,1 mV | A_VT ≈ 4,2–5 mV·µm → 18,0–21,7 mV | Pelgrom kuralı (Berkeley EE240B, özet) | ✗: 3,5–5,6 kat küçük | değişkenlik |
| Exp7 | ≥ 100 aygıt planı | — | — | yöntem | ölçüm planı |
| Exp8 | fisyon üssü −1,08 | damlalar Rayleigh sınırında bölünüyor | elektrodinamik terazi ölçümleri (özet) | ✓: %8 | bölünme ilkesi |
| Exp9 | araç çapraz denetimi | — | — | yöntem | doğrulama |
| Exp10 | 2B'de K < 0, bölünme yok | düşük yoğunlukta negatif sıkıştırılabilirlik | Eisenstein 1992, PRL 68, 674 (özet) | kısmen: var, eşik sayısı yok | kısıt |
| Exp11 | otokataliz gerekli | kendini çoğaltan noktalar deneyde | Lee 1994, Nature 369, 215 (özet) | ✓ | ilke |
| Exp12 | 4 K yüzer gövde kilitlenmesi, kink + histerezis | 180 nm NMOS'ta yalnız artan V_d'de kink, 30 K'de var, 40 K'de yok, PMOS'ta yok | Huang 2020, arXiv 2003.00207 (tam metin) | ✓ (V_M = 0,80 V ölçülmedi) | hücre kenedi |
| Exp13 | balast λ genişliği seçer | — | — | ölçüm yok | hücre genişliği |
| Exp14 | bölünme boyut penceresi | — | — | ölçüm yok | 1B şerit kısıtı |
| Exp15 | 1B delikler birleşir, B korur | manyetize Q-makinesinde elektron delikleri ve birleşme | Saeki 1979, PRL 42, 501 (atıflardan) | ✓ | iyon bulutu |
| Exp16 | PySINDy | — | — | yöntem | yöntem |
| Exp17 | j_on, w_s, N ≈ I/(0,75 j_on w_s) | Si'de yok; benzer: FEEP'te yayan iğne sayısı akımla artıyor | Mühlich 2022 (tam metin) | ölçüm yok (Si) | I_c, N(t) |
| Exp18 | +30 mV tohum 10/10 | tohum, ölçek eşleşmezliğinin 7,3–8,9σ'sı (10 µm × 0,15 µm) | Pelgrom ölçeği | tutarlı (dolaylı) | doğum yeri |
| Exp19 | V_DD adımı k = 1, 2, 4 | — | — | ölçüm yok | itki basamakları |
| Exp20 | zorla düşüş 3 → 7 | — | — | ölçüm yok | hızlı artış |
| Exp21 | 300 K'de bölünme yok | kink 40 K'de kayboluyor; n-GaAs filamanları 4,2 K'de | Huang 2020; Brandl 1989 (tam metin) | ✓ | kısıt: ≲ 30 K |
| Exp22 | 4 K'de 9,5e12'ye kadar tek alt bant; kaçak fA–pA | Si(100)'de ikinci alt bant ~6e12 üstünde; kaçak sayısı bulunamadı | Ray–Martin, IBM (özet) | ✗ alt bant (önceden kaydedilmemiş ek) / kaçak ölçüm yok | kayıp |
| Exp23 | aktarım 0,985–1,023 | — | — | ölçüm yok | η_aktarım |
| Exp24 | ~1,2 G₀ tek atomluk boğaz; ³¹P RF hedefi 1,000000 | HfO₂'de G₀ basamakları; ³¹P kapı > %98, okuma > %99,8, T2 60 ms | Long 2013, APL 102, 183505; Pla 2013, Nature 496, 334 (özet) | kısmen: var ✓, sayı iyimser | programlanabilir atom |
| Exp25 | GHZ 0,987 | Si GHZ %92,5 / %96,2 | Mądzik 2022; Thorvaldson 2025 (tam metin) | ✗: 2,5–6,2 puan iyimser | kuantum denetim |
| Exp26 | 17,5 s'de 0,86 | MREV-16 T2 = 25 s → 0,50 | Ladd 2005 (tam metin) | ✗: iyimser; "depolama çöker" sonucu daha da güçlü | kısıt: yükle denetim |
| Exp27 | F/P ≤ 1/c (yansımada 2/c); Casimir; DCE | IKAROS 1,12 mN ≤ 2SA/c = 1,78 mN (1 AU); Casimir %5; DCE fotonları | JAXA (özet); Lamoreaux 1997 (özet); Wilson 2011 (tam metin) | ✓ | karşılaştırma tabanı |

Exp28 denklemi Exp29'da NEXT'in 40 ölçülmüş noktasında doğrulandı (kök ortalama kare %1,45).

## Ölçülmüş atom yayımı itkisi (indiyum FEEP)

**Kaynak.** Mühlich vd., Acta Astronautica 197 (2022) 107–114, tez 3. bölüm (tam metin okundu). Aynı IFM Nano itkicinde:
- doğrudan itki: mN tork sarkacı, doğruluk < %2 (< 20 µN);
- dolaylı itki: T = I_em · √(2m·η_ν·V_em/e) · α · γ.

Bu denklem Exp28'in F = η_b · I · √(2mV/q) denklemiyle aynı.

**Ölçülmüş girdiler.**
- RPA ile iyon enerjisi tam olarak V_em (η_ν = 1).
- Akımın %98'i tek yüklü In⁺.
- γ, 160°'lik Faraday kabı taramasından.

**Sonuç (doğrudan − dolaylı):**

| koşul | fark |
|---|---|
| V_em = 5 kV sabit, nominal 2–4 mA | ortalama 4,4 µN (%1,6) |
| 3,5 mA'de 3–10 kV taraması | 8,5 µN (%2,8) |
| V_ex sabit | %4,3 (akım dalgalanmaları) |

İtki, denklemin öngördüğü gibi √V_em ve I_em ile ölçekleniyor. Önceden kayıtlı ölçüt (≤ %5): **✓**.

**"Programlanabilir" yön.** Yayım alanla açılıyor; iğneler belirli bir alan şiddetinde ateşleniyor ve yayan iğne sayısı akımla artıyor:
6 kV'ta 37 µA'de 28 iğneden 1'i, 3,5 mA'de 28'i. Exp18–19'un "hangi hücrenin, kaç hücrenin ateşleneceği" sorusunun ölçülmüş
karşılığı bu; silisyum hücrelerin kendisi değil.

**İndiyum buharlaşma alanı.** Exp28'in imaj tümseği değeri Q₀²/1,44 = 12,2 V/nm. Ölçülen LMIS uç alanı > 10 V/nm ile tutarlı; ölçülmüş
tek bir In değeri bulunamadı: **kısmen**.

## Ölçüme göre düzeltilmiş yol ve denklem

  F = γ · N(t) · I_c · √(2 m V / q),   yayım: alanla açılan atom yayıcılar,   tek yük: I = e·f

| halka | ölçülmüş dayanak | durum |
|---|---|---|
| itki = akım × √(2mV/q) × γ | NEXT 40 nokta (%1,45); In FEEP (%1,6–4,3) | ölçümle doğrulandı |
| yayımı alanla açmak, sayısı akımla artan yayıcılar | FEEP iğneleri (1/28 → 28/28); SiOx programlanabilir yayıcı (Exp29) | ölçüldü (benzer) |
| tek yük metrelemesi | Si pompası 0,27 ppm, T_e 7–9 K; istek üzerine tek iyon 123/139 | ölçüldü |
| turnike adası | ölçüm: ~2 nm'de 300 K, 5 nm'de < 200 K; model 3–4 kat kötümser | düzeltildi |
| hücre mekanizması sınıfı (darbe iyonlaşması + donmuş gövde) | 180 nm CMOS'ta 30 K'de var, 40 K'de yok | ölçüldü (sınıf), ≲ 30 K |
| uzayda kriyojenik çalışma | JWST MIRI kriyosoğutucuyla < 7 K (ESA) | ölçüldü |
| doğum yeri seçimi | +30 mV tohum ≈ 7–9σ ölçek eşleşmezliği | tutarlı (dolaylı) |
| hücre yasaları (λ, w_s, I_c, N basamakları, aktarım) | yok | **ölçülmedi** |

## Kalan ölçüm (laboratuvar; burada yapılamaz)

İtki, akım ve gerilimden ölçülmüş doğrulukla çıktığı için kalan katman yalnız elektriksel ölçüm istiyor. Önerilen sıra
(Exp12/13/17–20 önerileri):
1. **Geniş SKY130 nFET, ≤ 30 K** (tercihen 4 K) prob istasyonu. V_gs ≈ V_T − 0,1 V, seri dirençli drenaj; sabit V_ds bölümü
   (~0,80 V) ve yalnız artan taramada histerezis aranır. Huang 2020 ile aynı ölçüm türü.
2. **Bölümlenmiş drenaj + poly balast.** Bölüm akımları ayrı ölçülür; hücre başına akım I_c ve sayı basamakları
   N ≈ I/(0,75 j_on w_s) aranır.
3. **+30 mV tohumla doğum yeri** (Exp18); V_DD adımı ve zorla düşüşle k = 1, 2, 4 ve 3 → 7 (Exp19–20).
4. **Hücre dizisini bir In-FEEP ya da iyonik sıvı yayıcıya bağlamak.** İtki, doğrulanmış denklemle J ve V'den okunur.

## Kaynaklar

- Mühlich vd., Acta Astronautica 197 (2022) 107; tez: https://repositum.tuwien.at/handle/20.500.12708/136207 (tam metin)
- Shin vd., arXiv 1001.3724 (tam metin)
- Huang vd., arXiv 2003.00207 (tam metin)
- Zhao vd., arXiv 1703.04795 (tam metin)
- Shi vd., arXiv 1208.0519 (özet)
- Lenzlinger ve Snow, J. Appl. Phys. 40, 278 (1969) (özet)
- Pelgrom kuralı, Berkeley EE240B ders notu (özet)
- Eisenstein, Pfeiffer, West, PRL 68, 674 (1992) (özet)
- Lee, McCormick, Pearson, Swinney, Nature 369, 215 (1994) (özet)
- Saeki vd., PRL 42, 501 (1979) (atıflardan)
- Ray ve Martin, IBM (özet)
- Long vd., APL 102, 183505 (2013) (özet)
- Pla vd., Nature 496, 334 (2013); arXiv 1302.0047 (özet)
- IKAROS, JAXA 2010 (özet)
- Lamoreaux, PRL 78 (sayı 1), 1997 (özet)
- ESA, MIRI < 7 K (özet)
- Gabovich, Usp. Fiz. Nauk 140, 137 (1983) (tam metin)
- Önceki doğrulanmış kaynaklar: Ladd 2005, Mądzik 2022, Thorvaldson 2025, Wilson 2011, Brandl 1989 (tam metin)
