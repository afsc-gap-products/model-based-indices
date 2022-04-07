create or replace PROCEDURE UPDATE_SURVEY AS 
    CURSOR c_update IS
        (SELECT cruise_id,cruise, vessel_id, vessel_name,cruisejoin
         FROM race_data.v_cruises
         WHERE survey_definition_id IN (98, 47, 52, 143, 78, 112)
           AND year >= 1982
           AND cruisejoin IS NOT NULL
         MINUS
         SELECT cruise_id, cruise, vessel, UPPER(vessel_name), cruisejoin
         FROM survey);
    CURSOR c_haul_update IS
        (SELECT DISTINCT cruisejoin
        FROM SAFE.survey);
         
    v_cruisecount NUMBER;
    v_survey VARCHAR2(9);
    v_survey_id NUMBER;
    
BEGIN
    FOR r_update IN c_update
    LOOP
        SELECT COUNT(*)
        INTO v_cruisecount
        FROM survey
        WHERE cruise_id = r_update.cruise_id;
 
        SELECT DECODE(survey_definition_id,98,'EBS_SHELF',47,'GOA',52,'AI',143, 'NBS_SHELF',112, 'NBS_SHELF',78,'EBS_SLOPE',NULL)
        INTO v_survey
        FROM race_data.v_cruises
        WHERE cruise_id = r_update.cruise_id
        AND year >= 1982;   
        
        IF v_cruisecount = 0
        THEN
          INSERT INTO survey (survey,year,cruise,vessel, vessel_name,cruisejoin,cruise_id)
          VALUES (v_survey, SUBSTR(r_update.cruise,1,4), r_update.cruise, r_update.vessel_id,INITCAP(r_update.vessel_name), r_update.cruisejoin,r_update.cruise_id);
        ELSE
          UPDATE survey
          SET cruisejoin = r_update.cruisejoin,
              vessel = r_update.vessel_id,
              vessel_name = INITCAP(r_update.vessel_name)
          WHERE cruise_id = r_update.cruise_id;
          
        END IF;
    END LOOP;
    COMMIT;

    FOR r_haul_update IN c_haul_update
    LOOP
        -- Update EBS or NBS valid hauls table
        SELECT survey_definition_id
        INTO v_survey_id
        FROM race_data.v_cruises
        WHERE cruisejoin = r_haul_update.cruisejoin;
        
        IF v_survey_id = 98
        THEN
            MERGE INTO SAFE.ebs_vast_hauls s
            USING (SELECT * FROM racebase.haul WHERE cruisejoin=r_haul_update.cruisejoin AND haul_type = 3 AND performance >= 0 AND abundance_haul = 'Y') r
            ON (s.hauljoin = r.hauljoin)
            WHEN MATCHED THEN UPDATE SET
                s.CRUISE = r.CRUISE,
                s.HAUL = r.HAUL,
                s.HAUL_TYPE = r.HAUL_TYPE,
                s.PERFORMANCE = r.PERFORMANCE,
                s.START_TIME = r.START_TIME,
                s.DURATION = r.DURATION,
                s.DISTANCE_FISHED = r.DISTANCE_FISHED,
                s.NET_WIDTH = r.NET_WIDTH,
                s.NET_MEASURED = r.NET_MEASURED,
                s.NET_HEIGHT = r.NET_HEIGHT,
                s.STRATUM = r.STRATUM,
                s.START_LATITUDE = r.START_LATITUDE,
                s.END_LATITUDE = r.END_LATITUDE,
                s.START_LONGITUDE = r.START_LONGITUDE,
                s.END_LONGITUDE = r.END_LONGITUDE,
                s.STATIONID = r.STATIONID,
                s.GEAR_DEPTH = r.GEAR_DEPTH,
                s.BOTTOM_DEPTH = r.BOTTOM_DEPTH,
                s.BOTTOM_TYPE = r.BOTTOM_TYPE,
                s.SURFACE_TEMPERATURE = r.SURFACE_TEMPERATURE,
                s.GEAR_TEMPERATURE = r.GEAR_TEMPERATURE,
                s.WIRE_LENGTH = r.WIRE_LENGTH,
                s.GEAR = r.GEAR,
                s.ACCESSORIES = r.ACCESSORIES,
                s.SUBSAMPLE = r.SUBSAMPLE,
                s.ABUNDANCE_HAUL = r.ABUNDANCE_HAUL,
                s.AUDITJOIN = r.AUDITJOIN,
                s.CRUISEJOIN = r.CRUISEJOIN,
                s.REGION = r.REGION,
                s.VESSEL = r.VESSEL
            WHEN NOT MATCHED THEN 
            INSERT (CRUISE, 
                HAUL, 
                HAUL_TYPE, 
                PERFORMANCE, 
                START_TIME, 
                DURATION, 
                DISTANCE_FISHED, 
                NET_WIDTH, 
                NET_MEASURED, 
                NET_HEIGHT, 
                STRATUM, 
                START_LATITUDE, 
                END_LATITUDE, 
                START_LONGITUDE, 
                END_LONGITUDE, 
                STATIONID, 
                GEAR_DEPTH, 
                BOTTOM_DEPTH, 
                BOTTOM_TYPE, 
                SURFACE_TEMPERATURE, 
                GEAR_TEMPERATURE, 
                WIRE_LENGTH, 
                GEAR, 
                ACCESSORIES, 
                SUBSAMPLE, 
                ABUNDANCE_HAUL, 
                AUDITJOIN, 
                CRUISEJOIN, 
                HAULJOIN, 
                REGION, 
                VESSEL)
            VALUES (
                r.CRUISE,
                r.HAUL,
                r.HAUL_TYPE,
                r.PERFORMANCE,
                TO_DATE(r.START_TIME),
                r.DURATION,
                r.DISTANCE_FISHED,
                r.NET_WIDTH,
                r.NET_MEASURED,
                r.NET_HEIGHT,
                r.STRATUM,
                r.START_LATITUDE,
                r.END_LATITUDE,
                r.START_LONGITUDE,
                r.END_LONGITUDE,
                r.STATIONID,
                r.GEAR_DEPTH,
                r.BOTTOM_DEPTH,
                r.BOTTOM_TYPE,
                r.SURFACE_TEMPERATURE,
                r.GEAR_TEMPERATURE,
                r.WIRE_LENGTH,
                r.GEAR,
                r.ACCESSORIES,
                r.SUBSAMPLE,
                r.ABUNDANCE_HAUL,
                r.AUDITJOIN,
                r.CRUISEJOIN,
                r.HAULJOIN,
                r.REGION,
                r.VESSEL
                );
        ELSIF v_survey_id IN (143, 112)
        THEN
            MERGE INTO SAFE.nbs_vast_hauls s
            USING (SELECT * FROM racebase.haul WHERE cruisejoin=r_haul_update.cruisejoin AND haul_type = 3 AND performance >= 0 AND abundance_haul = 'Y') r
            ON (s.hauljoin = r.hauljoin)
            WHEN MATCHED THEN UPDATE SET
                s.CRUISE = r.CRUISE,
                s.HAUL = r.HAUL,
                s.HAUL_TYPE = r.HAUL_TYPE,
                s.PERFORMANCE = r.PERFORMANCE,
                s.START_TIME = r.START_TIME,
                s.DURATION = r.DURATION,
                s.DISTANCE_FISHED = r.DISTANCE_FISHED,
                s.NET_WIDTH = r.NET_WIDTH,
                s.NET_MEASURED = r.NET_MEASURED,
                s.NET_HEIGHT = r.NET_HEIGHT,
                s.STRATUM = r.STRATUM,
                s.START_LATITUDE = r.START_LATITUDE,
                s.END_LATITUDE = r.END_LATITUDE,
                s.START_LONGITUDE = r.START_LONGITUDE,
                s.END_LONGITUDE = r.END_LONGITUDE,
                s.STATIONID = r.STATIONID,
                s.GEAR_DEPTH = r.GEAR_DEPTH,
                s.BOTTOM_DEPTH = r.BOTTOM_DEPTH,
                s.BOTTOM_TYPE = r.BOTTOM_TYPE,
                s.SURFACE_TEMPERATURE = r.SURFACE_TEMPERATURE,
                s.GEAR_TEMPERATURE = r.GEAR_TEMPERATURE,
                s.WIRE_LENGTH = r.WIRE_LENGTH,
                s.GEAR = r.GEAR,
                s.ACCESSORIES = r.ACCESSORIES,
                s.SUBSAMPLE = r.SUBSAMPLE,
                s.ABUNDANCE_HAUL = r.ABUNDANCE_HAUL,
                s.AUDITJOIN = r.AUDITJOIN,
                s.CRUISEJOIN = r.CRUISEJOIN,
                s.REGION = r.REGION,
                s.VESSEL = r.VESSEL
            WHEN NOT MATCHED THEN 
            INSERT (CRUISE, 
                HAUL, 
                HAUL_TYPE, 
                PERFORMANCE, 
                START_TIME, 
                DURATION, 
                DISTANCE_FISHED, 
                NET_WIDTH, 
                NET_MEASURED, 
                NET_HEIGHT, 
                STRATUM, 
                START_LATITUDE, 
                END_LATITUDE, 
                START_LONGITUDE, 
                END_LONGITUDE, 
                STATIONID, 
                GEAR_DEPTH, 
                BOTTOM_DEPTH, 
                BOTTOM_TYPE, 
                SURFACE_TEMPERATURE, 
                GEAR_TEMPERATURE, 
                WIRE_LENGTH, 
                GEAR, 
                ACCESSORIES, 
                SUBSAMPLE, 
                ABUNDANCE_HAUL, 
                AUDITJOIN, 
                CRUISEJOIN, 
                HAULJOIN, 
                REGION, 
                VESSEL) 
            VALUES (
                r.CRUISE,
                r.HAUL,
                r.HAUL_TYPE,
                r.PERFORMANCE,
                TO_DATE(r.START_TIME),
                r.DURATION,
                r.DISTANCE_FISHED,
                r.NET_WIDTH,
                r.NET_MEASURED,
                r.NET_HEIGHT,
                r.STRATUM,
                r.START_LATITUDE,
                r.END_LATITUDE,
                r.START_LONGITUDE,
                r.END_LONGITUDE,
                r.STATIONID,
                r.GEAR_DEPTH,
                r.BOTTOM_DEPTH,
                r.BOTTOM_TYPE,
                r.SURFACE_TEMPERATURE,
                r.GEAR_TEMPERATURE,
                r.WIRE_LENGTH,
                r.GEAR,
                r.ACCESSORIES,
                r.SUBSAMPLE,
                r.ABUNDANCE_HAUL,
                r.AUDITJOIN,
                r.CRUISEJOIN,
                r.HAULJOIN,
                r.REGION,
                r.VESSEL
                );
        END IF;

    END LOOP;
    COMMIT;
    
END UPDATE_SURVEY;
/

BEGIN
  update_survey;
END;
/