BEGIN {FS=",";OFS = ";";lineCount=0;fs=250.0;ADS1299_gain = 24.0;offset_negative=15625*ADS1299_gain;negative_add_count=-2;ADS1299_Vref=4.5; scale_fac_uVolts_per_count = ADS1299_Vref / ((2^23)-1) / ADS1299_gain  * 1000000;acc1pre = 0; acc2pre = 0; acc3pre = 0;
print("index;time_sec;ringpar;ch1;ch2;ch3;ch4;ch5;ch6;ch7;ch8");
}
{ 
	if (NF >8){
	
	$10 = ++lineCount;
	$11 = $10/fs;
	$1 = sprintf("%d", strtonum("0x" $1));

	if(substr($2,1,1) > 7){ $2 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $2)+negative_add_count))}else{$2 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $2))};
	if(substr($3,1,1) > 7){ $3 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $3)+negative_add_count))}else{$3 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $3))};
	if(substr($4,1,1) > 7){ $4 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $4)+negative_add_count))}else{$4 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $4))};
	if(substr($5,1,1) > 7){ $5 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $5)+negative_add_count))}else{$5 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $5))};
	if(substr($6,1,1) > 7){ $6 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $6)+negative_add_count))}else{$6 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $6))};
	if(substr($7,1,1) > 7){ $7 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $7)+negative_add_count))}else{$7 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $7))};
	if(substr($8,1,1) > 7){ $8 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $8)+negative_add_count))}else{$8 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $8))};
	if(substr($9,1,1) > 7){ $9 = sprintf("-%f",   offset_negative - scale_fac_uVolts_per_count*(strtonum("0x" $9)+negative_add_count))}else{$9 = sprintf("%f",  scale_fac_uVolts_per_count*strtonum("0x" $9))};
	#if (1) {if(substr($10,1,1) == "F"){ $10 = sprintf("-%f",   strtonum("0x" "FFFF" $10));acc1pre = -1*strtonum("0x" $10);}else{$10 = sprintf("%f",  strtonum("0x" "0000" $10));acc1pre = strtonum("0x" $10);};} else {$10 = sprintf("%f", acc1pre)};	
	#if ($11) {$11 = sprintf("%f",  "0x" $11);acc2pre = strtonum("0x" $10);} else {$11 = sprintf("%f", acc2pre)};
	#if ($12) {$12 = sprintf("%f",  "0x" $12);acc3pre = strtonum("0x" $10);} else {$12 = sprintf("%f", acc3pre)};
        print($10,$11,$1,$2,$3,$4,$5,$6,$7,$8,$9);
        }
	
}