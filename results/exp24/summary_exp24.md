# Exp24 — Doğal atomları hücrelerle programlamak (0/1)

Tarih: 2026-10-15. Kod: `experiments/exp24_atom_program/`
- `ap24_scale.py`: hücre → atom ölçekleme yasası ve sınaması;
- `ap24_atoms.py`: hücre gövdesindeki atomlar, Coulomb camı, şekil;
- `ap24_catalog.py`: teknik kataloğu, denklemler ve sayılar.

Tahminler koşulardan önce yazıldı: `predictions_scale.json`, `predictions_atoms.json`. Şekil: `fig_exp24.png`.

## Soru ve kurallar

Kullanıcının hedefi: atomu dışarıdan programlamak, ona 0 ya da 1 olmasını söylemek. Önerisi: "hücreleri programlayabiliyorsak
hücreleri oluşturan şeyi de programlarız; hücre atom kadar küçülebilir"; atomu hücrenin içine hapsedip orada yönetmek.

"Atom programlama" için dört ölçüt kullandım:
- **(S)** iki kararlı durum;
- **(C)** hücreyle dış kontrol;
- **(A)** tek atom adreslenebilirliği;
- **(R)** okunabilirlik.

Bir atomun üç ayrı "durumu" programlanabilir: konumu, yük durumu ve spini. Bunlar ayrı ayrı belirtildi; biri diğerinin yerine
sayılmadı.

## Doğru sonuç (kısa)

1. **Silikon hücrenin kendisi atom boyuna inemez.** Türetilen tabanlar:
   - donmuş gövdenin bağlanma uzunluğu ℓ_b ≈ 2,2 µm (bu modelde);
   - tükenme genişliği W_dep ≈ 83 nm;
   - yüksek durumun 10 yıllık kararlılığı ≥ 36 nm;
   - atom boyutunda (0,5 nm) bir kutuda ortalama 2·10⁻⁵ katkı atomu var, yani gövde yok.

   Sınamada programlanabilirlik ℓ_b altında kayboldu. "Programlanabilen hücre atom kadar küçülür" önermesi silikonun
   elektronik filamanı için yanlış. Oksitteki iyonik filaman hücresi için ise doğru: tek atom genişliğine iniyor ve
   genişliğini silikon hücre seçiyor (§3c).
2. **Ama hücrenin içine hapsolmuş doğal atomlar tek tek 0/1 yapılabiliyor (model düzeyinde, ≤ 1 K'de güvenle).**
   - Hücre gövdesi donmuş bor atomlarıyla dolu bir kutu. Her atom ya deliği bağlı (B⁰) ya iyonize (B⁻).
   - Atom atom Coulomb camı hesabında her atomun kendi anahtarlama gerilimi ("adresi") var. 100 nm'lik bir gövdede yaklaşık
     94 atom 0 → 1 dönüyor.
   - Sonlu sıcaklık sınaması (Metropolis) sonrası, hücrenin kendi kararlı gövde gerilimi penceresinde (0,805–1,006 V):
     - 1 K'de 100 nm gövde başına ~35 atom her iki yanda ≥ 8 kT/e paya sahip; %88 başarıyla ~31 atom hücrenin kendi
       kontrolüyle güvenle 0/1 yapılıyor;
     - 4 K'de bu pencerede güvenle programlanabilen atom yok (0); ek gövde ön gerilimiyle ~2 atom.
   - Atomların %15–30'u histerezisli (kalıcı bit adayı); bu sayı sıfır sıcaklık dinamiğinden geldiği için üst sınır.
   - Kayıtlı ısıl tahminler (ortalama doğruluk ≥ 0,9) tutmadı: 4 K'de 0,82, 1 K'de 0,88.
   - Adresleme boyutla değil gerilimle yapılıyor.
   - **4 K için bulunan yöntem (sınandı, §2c):** gövde ~60 nm, adresler 4 K'de dengede kalibre, pay her yanda ≥ 8 kT/e.
     Böylece 60 nm gövde başına ~3,3 atom hücrenin kendi penceresinde ~%99 doğrulukla 0/1 yapılıyor (kayıtlı SB1 ✓).
     40 nm'de adresler pencerenin üstüne kayıyor (SB2 ✗).
   - Programlanabilen en küçük birim tek bir atom: deliğini bağlayan bor akseptörü (Bohr yarıçapı ~1,4 nm), ~60 nm'lik bir
     hücre gövdesinin içinde. Akım filamanı hücresinin kendisi ise atom boyuna inmiyor.
3. **Atomun konumunu ve spinini programlamak da simüle edildi.** Hücre burada adresleyici ve sürücü:
   - Konum: oksitteki iyonik filaman hücresi tek atom genişliğine iniyor. Silikon hücrenin gövde gerilimi genişliği seçiyor:
     0,81–0,90 V'ta 1 atom, 0,95–1,00 V'ta 2 atom (KMC, 20/20). Her yazma 2–4 oksijen atomunun yerini değiştiriyor (§3c).
   - Spin: hücre seçimi + tek RF darbesiyle tek bir ³¹P çekirdeği, P = 1, çapraz konuşma ≤ 1,7·10⁻⁷ (§3b).
   - Ek yol, literatür: akım darbesiyle tek atomluk metal temas (T4).
4. **Hücrenin kendi enerjisiyle örgü atomunu doğrudan oynatmak mümkün değil.** Exp23'ün sayıları geçerli: taşıyıcı enerjisi
   V_M'ye kilitli (0,82 eV).

## 1. Hücre atom boyuna inebilir mi? (ölçekleme yasası)

Denklemler:
- Hücre boyu: w ≥ 2λ_c, λ_c = √(2ΔV_c/(r_D j_on)). Rayı sertleştirmek (r_D büyük) hücreyi küçültür: r_D = 10 → 1e7 Ω/µm için
  λ_c = 36 µm → 36 nm.
- Programlama (doğum en büyük sürüşlü birimde, Exp18). w genişliğinde birim başına σ_u = A_VT/√(w L) (A_VT = 3,356 mV·µm,
  L = 0,15 µm):
  P_prog = P( max_{m tohum birimi}(a + σ_u z) > max_{N−m diğer}(σ_u z) ).
- Tabanlar:
  - W_dep = √(2ε·0,9 V/(q N_A)) = 83 nm;
  - gürültü (Kramers, doğum–ölüm yarı potansiyeli): ΔΦ = (c_b w/q)∫ln(i_çıkış/i_giriş)dV. 10 yıl için ΔΦ ≥ 26,5 gerekiyor;
    yüksek → düşük yönde bu w ≥ 36 nm demek (düşük → yüksek yönde 1,4 nm);
  - gövde bağlanması ℓ_b = √(L_gövde/(R_□ g)) ≈ 2,2 µm (R_□ = 1e12 Ω/□, g = 1e-13 S/µm). Bu taban sınamadan sonra fark edildi.

Sınama (`predictions_scale.json`; 5 durum × 10 alan; ilk doğumun tohumda olup olmadığı):

| λ_c | tohum | öngörü P_prog | kayıtlı aralık | gözlem | sonuç |
|---|---|---|---|---|---|
| 36 µm | 30 mV | 1,00 | 10–10 | 10/10 | ✓ |
| 3,6 µm | 30 mV | 0,86 | 7–10 | 10/10 | ✓ |
| 0,36 µm | 30 mV | 0,17 | 0–4 | 0/10 | ✓ |
| 0,036 µm | 30 mV | 0,06 | 0–2 | taramada doğum yok | değerlendirilemez |
| 0,36 µm | 150 mV | 0,997 | 10–10 | 0/10 | ✗ |

Son satırın nedeni sonradan bulundu: 2,2 µm'lik gövde bağlanması 0,36 µm'lik tohumu yayıp söndürüyor, ilk hücre akısız sınır
olan cihaz ucunda doğuyor. Kayıtlı formül "bağımsız birim" varsayıyordu ve λ_c < ℓ_b'de geçersiz.

Sonuç: hücre boyunun tabanı en az max(2λ_c, ℓ_b, W_dep, ~36 nm, litografi ~0,36 µm). Gerçek cihazda ℓ_b, kuyunun 4 K'deki
donma derecesine bağlı ve ölçülmeli.

**Mutlak taban (model parametrelerinden bağımsız).** Gövde tam donmuş (ℓ_b → 0) ve hücreler ayrık adalar (W_dep sınırı yok) olsa
bile iki temel sınır kalıyor:
- Isıl ömür (Kramers, yüksek → düşük bariyer integrali 0,117 V, c_b = 1 fF/µm): w ≥ ln(t/1 ms)·q/(c_b·0,117 V).
  1 ms için 3,2 nm; 1 s için 9,5 nm; 1 saat için 20,7 nm; 10 yıl için 36 nm.
- Kuantum hapsi: inversiyon kanalının yanal hapsetme enerjisi E₁ = ħ²π²/(2m*w²), m* = 0,19 m₀.
  - 10 nm: 19,8 meV;
  - 2 nm: 495 meV;
  - 1 nm: 1,98 eV, yani silisyum bant aralığının (1,17 eV) üstünde;
  - 0,5 nm: 7,9 eV.

  Bu yüzden ~1,3 nm'den dar bir elektronik kanal silisyumda var olamaz.

Kendiliğinden oluşan elektronik akım filamanının fiziksel alt sınırı ~1,3 nm, pratik sınırı ~3–10 nm (ms–s ömür). Atom ölçeği
(0,2–0,5 nm) dışarıda. Silisyumda atom genişliğinde iletken yollar yalnızca atomik dizilimle üretilmiş yapılarda var (örneğin
STM ile yazılmış Si:P teller); bunlar kendiliğinden oluşan, programlanabilen hücreler değil.

## 2. Hücrenin içine hapsolmuş atomlar (T1, simüle edildi)

**Model** (`ap24_atoms.py`):
- Rastgele konumlu bor atomları; N_A = 1,7e17 cm⁻³, ortalama aralık 18 nm.
- İnversiyon katmanı topraklı iletken düzlem.
- Büyük kanonik enerji:
  Ω = Σ n_i (E_a0 + eV_b + D_i) + ½ Σ_{i≠j} n_i n_j W_ij, burada:
  - n_i = 1 iyonize;
  - E_a0 = −E_g + E_A = −1,125 eV;
  - W_ij = k(1/r_ij − 1/r_ij*), k = e²/(4πε_Si) = 0,123 eV·nm;
  - D_i = −k/(4z_i) (öz-imge).
- Sıfır sıcaklıkta tek çevirme ve tek elektron sıçramasıyla gevşeme; V_b 0,25 mV adımla yukarı ve sonra aynı durumdan aşağı
  taranıyor.
- Kutular:
  - 100 nm ve 200 nm'lik yalıtılmış adalar (SOI benzeri, 150 nm kalın), her biri 5 tane;
  - 300 nm'lik periyodik "geniş cihaz", 2 tane.

4 K'de tükenme sınırı Debye uzunluğu kadar keskin (1,1 nm), atomlar ise 18 nm aralıklı. Bu yüzden sınır her atomu ayrı bir
gerilimde geçiyor.

**Sonuçlar:**

| kutu | atom | 0→0,97 V anahtarlanan | tek atomlu olay | tek başına adreslenebilir (≥ 5kT/e), 4 K / 1 K | histerezis > 5kT/e, 4 K / 1 K |
|---|---|---|---|---|---|
| 100 nm ada | ~245 | 94 | %64 | 43 / 96 atom | %15 / %30 |
| 200 nm ada | ~990 | 588 | %55 | 48 / 312 atom | %14 / %29 |
| 300 nm periyodik | ~2270 | 963 (taban yapaylığı dahil) | %42 | 1 / 347 atom | %17 / %27 |

**Tahminler (`predictions_atoms.json`):**
- A1 (eşik–derinlik sıra korelasyonu ≥ 0,8): yazıldığı haliyle ✗. İşareti yanlış yazmışım; fizikte derin atom önce nötrleşir,
  yani ρ < 0 olmalıydı. |ρ| ile:
  - periyodik cihazda (z < 130 nm) ✓, ρ = −0,95;
  - adalarda ✗, ρ = −0,35 ve −0,54; küçük adada eşik derinlikten çok yanal konuma da bağlı (üç boyutlu elektrostatik).
- A2 (1B adım formülü e·x_d/(εA), adalarda): ✗. 100 nm adada 1,8 mV ölçüldü, öngörü 5,3 mV. Küçük adalar V_b = 0'da
  tamamen tükenmiş, ilk nötrleşme ancak 0,69 V'ta başlıyor.
- A3 (anahtarlanan atom sayısı ±%15):
  - 100 nm ada ✓ (94'e karşı 99);
  - periyodik ✓ (963'e karşı 890; taban yapaylığı çıkarılınca ~850);
  - 200 nm ada ✗ (+%49).
- A4 (Coulomb etkileşimi merdiveni düzenler): ✗. Tersine, birden çok atomun birlikte döndüğü ilişkili olaylar tek başına
  adreslenebilen olayları azaltıyor.
- A5 (keşif): histerezis var; medyanı 100 nm adada 1,75 mV.

**Sınırlar:**
- Periyodik kutunun tabanında (z ≥ 130 nm) kutu başına ~113 sahte anahtarlama var; analizde ayrıldı.
- Sıfır sıcaklık ve tek hamle dinamiği histerezisi abartabilir.
- Dielektrik uyumsuzluğu ve S/D kenarları modelde yok.

**Tek atomu programlanabilir yapma ölçütü (türetilen):**
- k. atom, eşiği komşu eşiklerden δ ≥ m·kT/e (m ≈ 5) uzaksa 0/1 yapılabilir.
- Ortalama eşik aralığı ⟨ΔV⟩ ≈ e/C_etkin, C_etkin ≈ εA/x_d (1B).
- Poisson tabanında adreslenebilir kesir ≈ exp(−2mkT/(e⟨ΔV⟩)); gözlenen değer bunun biraz altında.
- 4 K'de iyi adresleme için hücre gövdesi ~(100 nm)² ölçeğinde olmalı ya da sıcaklık ≤ 1 K olmalı.

**Isıl sınama (`predictions_thermal.json`, `thermal.json`; koşudan önce kaydedildi):**
- Yöntem: 100 nm adalarda, 5 kutu × 8 atom. Her atom için gövde gerilimi eşiğin altındaki ve üstündeki boşluğun ortasına
  kondu; 4 K ve 1 K'de büyük kanonik Metropolis (tek çevirme + sıçrama, 1500 tarama) ile hedef durumun olasılığı ölçüldü.

| | 4 K | 1 K |
|---|---|---|
| ortalama P(hedef) (öngörü ≥ 0,9) | 0,82 ✗ | 0,88 ✗ |
| P ≥ 0,9 olan (atom, durum) | %56 | %69 |
| iki durumda da ≥ 0,9 olan atom | %28 | %47 |

- Doğruluk payla artıyor. 4 K'de yarım boşluk (pay):
  - 2,5–4 kT/e → 0,74;
  - 4–8 kT/e → 0,86;
  - 8–16 kT/e → 0,97 (durumların %88'i ≥ 0,9).
- 1 K'de büyük paylarda da arada başarısızlık var (8–16 kT/e → 0,81). Bu, sıfır sıcaklık taramasından okunan adresin
  dengedeki adresle her zaman örtüşmemesinden geliyor (cam davranışı, çok-cisim yeniden düzenlenmesi).
- **Düzeltilmiş ölçüt:** her yanda pay ≥ ~8 kT/e, artı sonlu sıcaklıkta kalibre edilmiş adres (yaz–doğrula).
- Hücrenin kendi kararlı yüksek durum gövde gerilimi V_D ile 0,805–1,006 V arasında ayarlanıyor. Bu pencerede, ölçütü
  karşılayan atom sayısı 100 nm gövde başına:
  - 1 K'de ~35 (26–40);
  - 4 K'de 0 (pencere dışında ~2).
- 4 K için gerekli gövde: aralık ⟨ΔV⟩ ∝ x_d/(εA). ≥ 8 kT/e iki yanlı pay için ⟨ΔV⟩ ≳ 16 mV gerekiyor; bu da 1B ile A ≈ (57 nm)²,
  adalardaki 3B azalmayla ~(35 nm)² demek (öngörü, simüle edilmedi). Hücrenin 1 K'deki çalışması modellenmedi (4 K modeli
  kullanıldı).

## 2c. 4 K'de çözüm: küçük gövde + dengede kalibre adres (sınandı)

`ap24_small.py`, `small_bodies.json`; tahminler koşudan önce yazıldı (`predictions_small.json`).

**Yöntem:**
- 40 ve 60 nm'lik adalar (10'ar), aynı Coulomb camı modeli.
- V_b 4 K'de 0,5 mV adımla Metropolis ile taranıyor, durum adımdan adıma taşınıyor. Her atomun ⟨n⟩(V_b) eğrisinin 0,5'i
  kestiği yer dengedeki adresi veriyor.
- Komşu adreslere her iki yanda ≥ 8 kT/e uzaklıkta olan atomlar, alt ve üst boşluğun ortasında programlanıyor; hedef durum
  olasılığı ayrı bir Metropolis koşusuyla ölçülüyor.

**Sonuçlar:**

| gövde | atom | dengede adreslenen | nitelikli (≥ 8 kT/e) | hücre penceresinde (0,805–1,006 V) | ortalama P(hedef) | P ≥ 0,9 olan çift |
|---|---|---|---|---|---|---|
| 40 nm | 34,9 | 6,1 | 2,1 | 0 | 0,996 | %100 (42 çift) |
| 60 nm | 82,6 | 13,5 | 7,2 | 3,3 | 0,991 | %99,3 (144 çift) |

**Tahminler:**
- SB1 (ortalama ≥ 0,9 ve çiftlerin ≥ %80'i ≥ 0,9): ✓.
- SB2 (40 nm'de pencerede ≥ 2 atom): ✗, 0 çıktı. Küçük adada adresler 1,006 V'un üstüne kayıyor; hücrenin kendi penceresinin
  dışında kalıyor. Bunlar için ek gövde ön gerilimi gerekir.
- SB3 (40 nm'de nitelikli oran daha yüksek): ✗. %6,0'a karşı 60 nm'de %8,7. Küçük adada dengede tek kesişimli atomlar az.

**Sonuç:** 4 K'de hücrenin kendi kontrolüyle tek atom 0/1 için doğru yöntem ~60 nm gövde, dengede kalibre adres ve ≥ 8 kT/e
pay. Böylece gövde başına ~3 atom ~%99 doğrulukla programlanıyor. 100 nm'deki kalibresiz %82 bu yöntemle %99'a çıktı. Önceki
TH1 tahmini (✗) kayıtta olduğu gibi kalıyor; bu yeni ve ayrı bir sınama.

**Okuma (tahmin):**
- Bir atomun yük değişimi kanalda yerel bir eşik kayması yapar: ΔV_T ≈ e/(C_ox A) × geometri ≈ 0,6 mV (100 nm, z ≈ 40 nm).
- 4 K'de eşik altı eğimi 20 mV/dek olduğundan bu ~%7 akım basamağı demek. Ölçülebilir; deneyde doğrulanmalı.

## 3. Teknik kataloğu (denklem + sayı, `catalog.json`)

| teknik | denklem | bizim sayılar | S/C/A/R |
|---|---|---|---|
| T1 hücre gövdesinde tek akseptör yükü | Ω(n) yukarıda | 100 nm gövdede 4 K'de ~43, 1 K'de ~96 adreslenebilir atom | ✓ ✓ ~ ~ |
| T7 hücre kenarından tünelleme kabuğu | Γ(d) = Γ₀ e^{−2d/a*}, a* = 1,37 nm | 1e3 /s'den 1e-3 /s'ye kabuk 14 → 24 nm, kalınlık 9,4 nm < atom aralığı 18 nm, yani bir anda ~1 atom katmanı | ✓ ✓ ~ ~ |
| T3 tek ³¹P çekirdek spini (n-kuyulu hücre) | ν(D⁺) = γ_n B; ν(D⁰) = A/2 ∓ γ_n B | B = 0,5–1,5 T'de ayrım 7–42 MHz ≫ kHz Rabi; tek RF tonu yalnız hücrenin seçtiği yük durumundaki atomu döndürür | ✓ ✓ ✓* ✓* |
| T4 tek atomluk metal temas (Schirm 2013) | \|I_kol\| = w·j_on ≥ I_s | I_s 10–50 µA → 3,9–19,5 µm'lik hücre kolu; iki yön için iki hücre | ✓ ✓ ✓ ✓ |
| T5 SKY130 ReRAM O-boşluğu zinciri | n = I_C/(V_C G₀); I_C yazma FET'inin gövde gerilimindeki akımı | W 1 µm: doğruluk %99,2, n ≈ 5,6; W 0,42 µm, V_T 0,70: %98,3, n ≈ 1,22 (p10–p90 1,06–1,46); tekme 0,20–0,26 V < pay 0,42 V | ✓ ✓ ✓ ✓ |
| T8 akseptör spini (EDSR + hücre alanı) | f_L = gμ_B B/h; Δf = δg μ_B B/h | δg = 1e-3–1e-2 → 4–140 MHz adres kayması; Rabi 3–40 MHz | ✓ ~ ~ ✓* |
| T2 cephe alanıyla iyonlaştırma | F_c = πεE_b²/e³; Γ(F) = (4/τ_a)(F_a/F)e^{−2F_a/(3F)} | 1 /s için 11,6 kV/cm; 0,5 µm cephe 18 kV/cm (hızlı), 2,2 µm cephe 4 kV/cm (yok) | ✓ ~ ✗ ~ |
| T6 Si–H bağı (titreşim merdiveni) | oran ~ (Φ p τ_v)^N / τ_v, N = 14 | Φ = 7e9 /s; günde 1 kopma için p_in ≥ 2e-3–1e-2 (bilinmiyor); hangi H'nin koptuğu rastgele | ✓ ✗ ✗ ~ |
| hücrenin kendisi atom boyunda | tabanlar §1 | ℓ_b, W_dep, gürültü; atom kutusunda 2e-5 katkı | ✗ |
| örgü atomunu doğrudan oynatma | Exp23 eşikleri | 2–11 mertebe eksik | ✗ |

(✓* = ek RF ve tek elektron okuma altyapısıyla.)

Literatür (bu oturumda doğrulanan):
- SkyWater SKY130 ReRAM: `sky130_fd_pr_reram`, metal1–metal2 MIM. Şekillendirme 2,2–3,1 V. SET: WL 1,7 V, BL 2,4 V.
  RESET: WL 2,5 V, SL 2,6 V. Darbeler 1 µs.
- HfO₂ ReRAM 4 K'de: SET 1,30 → 1,46 V, RESET 1,75 → 1,86 V, 1M döngü (Lu vd., Adv. Electron. Mater. 2026). 1,45 K'de
  şekillendirme ve anahtarlama (Pérez-Bosch Quesada vd., IEEE EDL 45, 2391, 2024).
- Kuantize iletkenlik ve tek atomluk zincir: Long vd., Appl. Phys. Lett. 102, 183505 (2013); 16 kuantize seviye: Chen vd.,
  ACS Appl. Electron. Mater. (2019).
- Akımla tek atom belleği: Schirm vd., Nature Nanotechnol. 8, 645 (2013). Anahtarlama 10–50 µA ve 50–200 mV, 250 mK'de
  >500 döngü.
- Ticari CMOS transistörde iki doğal bor atomunun spin okuması: van der Heijden vd. (2018). Si boşluk spininde EDSR,
  0,41 mV ile 3–40 MHz: Crippa vd. (2018).
- Bellekten, bu oturumda doğrulanmadı: ³¹P için A = 117,53 MHz ve γ_n = 17,235 MHz/T; iyonize ³¹P çekirdeğinin 28Si'de
  T2 = 35,6 s (Muhonen vd. 2014); RRAM'de R ≈ V_C/I_C, V_C ≈ 0,4 V (Ielmini 2011).

## 3b. Tek bir doğal atomun çekirdek spinini programlamak (T3, simüle edildi)

`ap24_spin.py`, `spin.json`; tahminler koşudan önce yazıldı (`predictions_spin.json`).
- Kurgu: hücre gerilim adresiyle hedef ³¹P katkı atomunu iyonize (D⁺), diğerlerini nötr (D⁰) bırakıyor. Bunun için n-gövdeli
  ayna hücre gerekiyor; varsayım, simüle edilmedi.
- Tek bir global RF darbesi f = γ_n B frekansında veriliyor. Hesap dönen dalga yaklaşımında iki seviyeli Rabi formülüyle ve
  sayısal tümlevle yapıldı.
- Sonuç:
  - B = 1 T, Rabi 10 kHz: π darbesi 50 µs; hedef çekirdek P = 1,000000 ile dönüyor;
  - nötr komşulardaki en büyük istenmeyen dönme 2,9·10⁻⁸ (m_s = −½) ve 1,7·10⁻⁷ (m_s = +½);
  - SP1 ✓ (≤ 10⁻⁶);
  - formül ile sayısal tümlev 10⁻¹⁵ düzeyinde örtüşüyor, SP2 ✓;
  - 1,5 T ve 100 kHz'de çapraz konuşma 2·10⁻⁴'e çıkıyor (tasarım seçimi).
- Anlamı: hücre doğru yük düzenini kurduğu sürece tek bir doğal atomun çekirdek spini neredeyse kusursuz biçimde 0 → 1
  yapılıyor. Toplam doğruluk ≈ (yük düzeni doğruluğu) × 1. Darboğaz spin değil, yük adreslemesi (§2).
- Okuma (Morello 2010 tipi spinden yüke dönüşüm) ve RF hattı literatür düzeyinde; burada simüle edilmedi.

## 3c. Tek atom genişliğinde programlanabilir hücre: oksitte iyonik filaman (KMC, simüle edildi)

`ap24_ionic.py`, `ionic.json`; tahminler koşudan önce yazıldı (`predictions_ionic.json`).

**Neden bu yol:** Silikondaki elektronik filaman ~0,1 µm'nin altına inemiyor (§1). Aynı tür hücre (bistabil ortamda yerel akım
kanalı, dışarıdan doğurulup söndürülebilen) oksitte iyonik olarak oluşuyor. Orada en küçük ölçek atomlar arası sıçrama
mesafesi, dolayısıyla hücre tek atom genişliğinde olabilir. SKY130'da bu katman (HfO₂ ReRAM) mevcut.

**Model:**
- 2B kafes kesiti: 9 × 16 site, a = 0,3 nm, 4,8 nm kalınlık. Her site ya O ya oksijen boşluğu (V_O).
- İletim: tabana bağlı boşluk kümesi tepeye ulaşırsa G = n G₀, n = en dar kesitteki boşluk sayısı (bir sütun = bir kanal;
  Long 2013). Ulaşmazsa d katmanlık boşluktan tünelleme.
- Olaylar: kümeye komşu bir O sitesi, ν·exp(−max(0, E_a − αF)/(kT_yerel)) hızıyla V_O'ya dönüşüyor.
  - F: uçtaki boşluk alanı, bağlı filamanda en dar kesit alanı.
  - T_yerel: Joule ısınmasıyla.
- Devre: bit hattı rampası → ReRAM → yazma nFET'i (4 K, EKV). Transistörün kapısı silikon hücrenin gövde gerilimi; akım
  sınırı I_C'yi bu belirliyor.

**Kalibrasyon:**
- α, doğrulanmış 4 K SET gerilimine (1,46 V, Lu 2026) göre ayarlandı.
- Varsayımlar: κ (tipik HfO₂ HRS ~100 kΩ, SkyWater değeri değil), ν = 1e13 /s, E_a = 1,0 eV, R_th = 1e6 K/W.
- Koşudan önce türetilen yasa: n = max(1, I_C/(G₀V_C)), V_C = (E_a − kT ln νt)·a/α ≈ 0,52 V.

**Sonuçlar** (4 K, 20'şer deneme; her koşuda SET 1,466 V):

| hücrenin gövde gerilimi | I_C | n (20 deneme) | öngörü | yer değiştiren O atomu |
|---|---|---|---|---|
| 0,81 V | 16,8 µA | 1 (20/20) | 1 | 2 |
| 0,85 V | 31,2 µA | 1 (20/20) | 1 | 2 |
| 0,90 V | 55,5 µA | 1 (20/20) | 1 | 2 |
| 0,95 V | 86,7 µA | 2 (20/20) | 2 | 4 |
| 1,00 V | 124,9 µA | 2 (20/20) | 3 | 4 |

**Tahminler:**
- IC1 (±1 içinde): ✓, 5 noktanın 4'ünde birebir.
- IC2 (I_C < G₀V_C = 40 µA'de tek atomluk boğaz ≥ %80): ✓, %100.
- IC3 (yer değiştiren atom ≤ boşluk + n + 2): ✓.
- IC4 (hücre penceresinde n 1'den ≥ 3'e): ✗. En fazla 2.

**Duyarlılık (V_b = 0,95 V):** E_a 0,8 ya da 1,2 eV → n = 2; R_th 1e5 → n = 1; R_th 1e7 → n = 3. Yani ısıl direnç ±1 kanal
belirsizlik getiriyor.

**Anlamı:**
- Programlanabilir bir akım filamanı hücresi tek atom genişliğine iniyor. Silikon hücre, kendi gövde gerilimiyle bu hücrenin
  1 mi 2 atom mu olacağını seçiyor. Yazma 2–4 oksijen atomunun konumunu değiştiriyor.
- Bu, "hücreleri programlayabiliyorsak onları oluşturan şeyi de programlarız" zincirinin doğru fiziksel karşılığı: silikon hücre
  (µm) → iyonik hücre (1–2 atom).

**Sınırlar:**
- 2B kafes soyutlaması.
- 4 K'de hızlar neredeyse eşik gibi davrandığı için sonuç deterministik; gerçek aygıtta istatistiksel dağılım beklenir (Long 2013
  histogramları).
- Varsayılan parametreler; ölçüm yok.

## 4. Sınırlar ve sonraki adım

- T1 bir model hesabı. Ölçüm önerisi:
  - (100 nm)² gövdeli SKY130 yapıda, 4 K ve 1 K'de gövde (ya da kanal) geriliminin ince taranması;
  - tek atom basamaklarının (~%7 akım) ve histerezisin aranması.
- T5 yarı statik hesaplandı; yazma darbesinin hücre dinamiğine etkisi zaman alanında simüle edilmedi, yalnız tekme denklemiyle
  hesaplandı.
- T3, T4 ve T8 için RF hattı, manyetik alan, tek elektron okuyucu (SET) ya da atomik temas üretimi gerekir. Burada denklem ve
  literatür düzeyinde değerlendirildi.
- Hücrenin kendisi atom boyuna inemez. Ama atom ölçeğindeki hedefler µm ölçeğindeki hücrelerle gerilim, frekans, akım ya da
  tünelleme mesafesi üzerinden tek tek adreslenebilir. Adresleme boyutla değil enerji ekseninde yapılıyor; MR cihazının spinleri
  frekansla seçmesine benzer.
